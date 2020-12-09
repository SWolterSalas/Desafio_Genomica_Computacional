##### Documentacion de script utilizado para el analisis metatranscriptomico de genomas de SARS-CoV-2

##### Autor: Sebastian Wolter Salas
##### <s.woltersalas@uandresbello.edu>
##### v 1.0
##### 20201208

##### Requisitos:
## Conda 

##### En un ambiente, instalar estos programas en el siguiente orden (aunque es recomendado instalarlas en ambientes separados):
##### Programas                Versiones            Instalar
#####
##    fastqc                   0.11.9               conda install -c bioconda fastqc
##    multiqc                  1.9                  conda install -c bioconda multiqc
##    trim-galore              0.6.6                conda install -c bioconda trim-galore
##    bowtie2                  2.4.2                conda install -c bioconda bowtie2
##    samtools                 1.11                 conda install -c bioconda samtools
##    bcftools                 1.8                  conda install -c bioconda bcftools
##    seqtk                    1.3                  conda install -c bioconda seqtk
##    mafft                    7.475                conda install -c bioconda mafft
##    iqtree                   1.6.12               conda install -c bioconda iqtree
##    pangolin                 1.1.14               conda install -c bioconda pangolin
##    llama                    0.1                  https://github.com/cov-lineages/llama
#####

##### Es necesario crear 2 directorios con los siguientes nombres:
## /1_fastas            Este directorio contendra todos los fastas para el analisis
## /1_referencia        Este directorio contendra el archivo de referencia. Por defecto, posee la referencia NC_045512.2.fasta, por lo que si es necesario cambiar, ver lineas 64, 93, 130, 145 y 150 para cambiarlo por aquella que necesites

##### Parte 1: Analisis de calidad

#### En primer lugar, se utiliza fastq para visualizar las estadisticas y calidad de las reads obtenidas mediante la secuenciacion de los metatranscriptomas

## Se utilizan 2 INPUTS, uno dedicado para los archivos fastq y otro para los archivos fq
INPUT1=$1 #"raw__raw_reads.read_2.fastq"
INPUT2=$2 #"raw__raw_reads.read_2.fq"

cd ..
mkdir 2_trimming  # Creacion de directorios de trabajo
mkdir 3_maping
mkdir 4_fastas_depurados
mkdir stats
mkdir 5_analisis_alineamiento
mkdir 6_llama
cd 1_fastas

fastqc ./* -t 16  # Analisis de los fastas contenidos en la carpeta 
multiqc ./
## Se identifican los problemas de calidad y se procede a realizar el trimming mediante Trim-Galore!. Debido al buen estado de las reads, es que se filtro por parametros relativamente exigentes

##### Parte 2: Trimming y analisis de calidad

## Se realiza el trimming mediante un ciclo para todos los archivos fasta de manera apareada
for archivo in $INPUT1
do
	TEXTO=${archivo%.raw__raw_reads.read_*}
	trim_galore -q 25 --o /home/swolter/Downloads/desafio_genomica/2_trimming --length 50 --trim-n -j 16 --paired $TEXTO\.raw__raw_reads.read_1.fastq $TEXTO\.raw__raw_reads.read_2.fastq
done

cd ../1_referencia

## Se construyen los indices de la referencia mediante Bowtie2
bowtie2-build ref_NC_045512.2.fasta referencia.index

cd ../2_trimming

## Se evalua la calidad nuevamente, con tal de determinar si el trimming permitio aumentar la calidad de las reads
fastqc ./* -t 16
multiqc ./

##### Parte 3: Mapeo

## Se realizo el mapping mediante Bowtie2, utilizando los pares de los fastas correspondientes
for archivo in $INPUT2
do
	TEXTO=${archivo%.raw__raw_reads.read_*}
	bowtie2 -p 16 -x ./../fastq/referencia/referencia.index -1 $TEXTO\.raw__raw_reads.read_val_1.fq -2 $TEXTO\.raw__raw_reads.read_val_2.fq -S ./../3maping/$TEXTO\.raw__raw_reads.read.sam | samtools view -bS | samtools sort > $TEXTO\.raw__raw_reads.read_val_2.bam
done 

##### Parte 4: Consenso de variantes y formacion de archivos fastas

## El mapeo anterior, tendra de output archivos bam, por lo que es necesario adaptar el nombre de los archivos de INPUT a bam
INPUT3=$3 #"BC*_val_2.bam"

cd ../3_mapping

## Para cada archivo generado mediante el mapeo, se utiliza bcftools para generar un consenso y llamar a las variantes mas preponderantes.
## Ademas, se utiliza vcfutils y seqtk para generar el archivo fasta consenso
for archivo in $INPUT3
do
	TEXTO=${archivo%.raw__raw_reads.read*}
	bcftools mpileup -f /home/swolter/Downloads/desafio_genomica/1_referencia/ref_NC_045512.2.fasta $TEXTO\.raw__raw_reads.read_val_2.bam > $TEXTO\.raw__raw_reads.read_2_val_2.vcf
	bcftools call -c $TEXTO\.raw__raw_reads.read_2_val_2.vcf | vcfutils.pl vcf2fq | seqtk seq -aQ64 -q20 > $TEXTO\.raw__raw_reads.read_val_2.fasta
	mv $TEXTO\.raw__raw_reads.read_val_2.fasta /home/swolter/Downloads/desafio_genomica/4_fastas_depurados/
done 

## Ahora, es necesario determinar las estadisticas, como por ejemplo el de la cobertura
## Es por esto que se convoca stats para cada archivo bam y se transforma en un txt, siendo movida a la carpeta stats
for archivo in $INPUT3
do
	TEXTO=${archivo%.raw__raw_reads.read*}
	samtools stats $TEXTO\.raw__raw_reads.read_val_2.bam > $TEXTO\.raw__raw_reads.read_val_2.txt
	mv /home/swolter/Downloads/desafio_genomica/3_mapping/$TEXTO\.raw__raw_reads.read_val_2.txt /home/swolter/Downloads/desafio_genomica/stats/$TEXTO\.raw__raw_reads.read_val_2.txt
done

## Creo el archivo que tendra todas las coberturas
echo ' ' > coverage.txt

## Para cada archivo bam leo la profundidad de las reads, utilizando estos parametros para determinar una cobertura aproximada
for archivo in $INPUT3
do
	TEXTO=${archivo%.raw__raw_reads.read*}
	echo "${TEXTO}" >> coverage.txt
	samtools depth -a $TEXTO\.raw__raw_reads.read_val_2.bam | awk '{sum+=$3} END { print "Average = ",sum/NR}' >> coverage.txt
done
mv /home/swolter/Downloads/desafio_genomica/3_mapping/coverage.txt /home/swolter/Downloads/desafio_genomica/stats/coverage.txt

##### Parte 5: Linajes, alineamientos y variantes

## El consenso generara un archivo fasta, por lo que es necesario adaptar el INPUT para la lectura de archivos fasta
INPUT4=$4 #"BC*.fasta"

cd ../4_fastas_depurados

## Para proseguir el con los analisis posteriores, es que se generan archivo fasta con el encabezado correcto, ya que en los pasos anteriores se modifica por el de la referencia
for archivo in $INPUT4
do
	TEXTO=${archivo%.raw__raw_reads.read*}
	cat $TEXTO\.raw__raw_reads.read_val_2.fasta | sed 's/\NC_045512.2/\'"$TEXTO\.raw__raw_reads.read_val_2.fasta"'/' > $TEXTO\\.raw__raw_reads.read_val_2.fasta2
done

## Se crea el archivo fasta_compilado, el cual contendra todos las secuencias de los analisis anteriores, en conjunto a su encabezado correcto
echo " " > fasta_compilado.fasta

## Se lee cada archivo y se agrega al documento creado previamente
for archivo in $INPUT4
do
	TEXTO=${archivo%.raw__raw_reads.read*}
	cat $TEXTO\\.raw__raw_reads.read_val_2.fasta2 >> fasta_compilado.fasta
done 

## Se mueve todo a otra carpeta de trabajo, copiando la referencia
mv fasta_compilado.fasta /home/swolter/Downloads/desafio_genomica/5_analisis_alineamiento/
cp /home/swolter/Downloads/desafio_genomica/1_referencia/ref_NC_045512.2.fasta /home/swolter/Downloads/desafio_genomica/5_analisis_alineamiento/

cd ../5_analisis_alineamiento

## Se realiza el alineamiento con MAFFT, el arbol filogenetico con IQTree y el analisis de variantes con Pangolin mediante el compilado creado en el paso anterior
mafft --thread 16 --auto --keeplength --addfragments fasta_compilado.fasta ref_NC_045512.2.fasta > alignment.fasta
pangolin fasta_compilado.fasta --outfile fasta_compilado.csv --write-tree -t 16
iqtree -s alignment.fasta

cd ../6_llama
# Muevo todos los archivos a una carpeta de trabajo nueva
mv /home/swolter/Downloads/desafio_genomica/5_analisis_alineamiento/alignment.fasta /home/swolter/Downloads/desafio_genomica/6_llama
mv /home/swolter/Downloads/desafio_genomica/5_analisis_alineamiento/fasta_compilado.csv /home/swolter/Downloads/desafio_genomica/6_llama
mv /home/swolter/Downloads/desafio_genomica/5_analisis_alineamiento/alignment.fasta.treefile /home/swolter/Downloads/desafio_genomica/6_llama

## Cambio los nombres de los archivos de entrada a llama, ya que no reconoce otros archivos que no posean esos nombres
mv alignment.fasta.treefile global.tree
mv fasta_compilado.csv metadata.csv

## Se realiza el analisis de SNPs desde el alineamiento hecho previamente mediante MAFFT, utilizando el arbol filogenetico de IQTree y el analisis de las cepas de Pangolin
llama -i metadata.csv -f alignment.fasta -d /mnt/c/Users/Wolter/Downloads/desafio_genomica/ --data-column name --lineage-representatives 
