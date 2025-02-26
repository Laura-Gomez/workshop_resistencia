# Taller "Genómica aplicada para fortalecer la lucha contra las bacterias multirresistentes"

## Objetivo

En este taller conocerás herramientas bioinformáticas para la identificación de genes de resistencia a antibióticos.

Para facilitar su implementación como parte de un flujo de trabajo estandarizado y reproducible utilizaremos Nextflow y Docker.

## Conectándonos al servidor y preparando el directorio de trabajo

**Para usuarios MAC:** Abre tu Terminal

**Para usuarios Windows:** Abre el Subsistima de Linux en Windows o abre Mobaxterm


Para conectarnos al servidor podemos usar el siguiente comando:

```
ssh alumnoN@10.0.15.11
```

Después de ingresar el comando anterior, teclea tu contraseña. **IMPORTANTE**: No verás nada en pantalla mientras ingresas la contraseña, aún así la contraseña se está escribiendo.


Una vez dentro del servidor, ubica tu directorio de trabajo:

```
pwd
```

Clona el repositorio con el código para este taller:

```
git clone https://github.com/Laura-Gomez/workshop_resistencia.git
```


## Aprendiendo a usar Nextflow

Muévete al directorio *workshop_resistencia/exercise_fastp*

```
cd workshop_resistencia/exercise_fastp
```

En este ejercicio ejecutaremos el comando fastp. Visualiza la ayuda de este programa 

```
fastp -h
```

Los datos de secuenciación se encuentran en el directorio */home/lgomez/workshop_resistencia/data*

```
ls /home/lgomez/workshop_resistencia/data
```

Tomando en cuenta que tenemos datos de secuenciación pareados, el comando fastp se ejecutaría de la siguiente manera para la muestra XX:

```
mkdir fastp_alone
 
fastp \
        --in1 /home/lgomez/workshop_resistencia/data/YV01_S2_R1_001.fastq.gz \
        --in2 /home/lgomez/workshop_resistencia/data/YV01_S2_R2_001.fastq.gz \
        --out1 fastp_alone/YV01_R1.fastq.gz \
        --out2 fastp_alone/YV01_R2.fastq.gz \
        --html fastp_alone/YV01.html \
        --json fastp_alone/YV01.json
```

Ahora vamos a trasladar este comando a un flujo de trabajo de Nextflow, para esto necesitamos:

**1. Archivo main.nf**
Contiene los pasos del flujo de trabajo, el orden de los módulos que se van a ejecutar

**2. Archivo modules.nf**
Contiene las instrucciones específicas para cada módulo del flujo de trabajo

**3. Archivo nextflow.config**
Contiene los parámetros del flujo de trabajo

Visualiza el contenido de estos archivos:
```
more main.nf
more modules.nf
more nextflow.nf
```

Ahora, ejecutemos el flujo de trabajo:

```
nextflow run main.nf
```

## Aprendiendo a conectar Nextflow y Docker

Muévete al directorio *workshop_resistencia/exercise_fastp*

```
cd ../workshop_resistencia/exercise_fastp
```

En esta parte de la práctica usaremos el comando bwa para hacer un alineamiento con las lecturas procesadas

Ve la ayuda del comando BWA

```
bwa -h
```

Este comando no está instalado en el servidor. Para ejecutarlo tenemos varias opciones, dos de las cuales son:

1. Instalarlo [https://bio-bwa.sourceforge.net/]
2. Usarlo a través de un contenedor de Docker

¡Nextflow está preparado para recibir las indicaciones de cual contenedor usar para cada uno de los modulos!. Solo hay que descargar el Docker que nos interesa, agregar la información del contenedor a usar junto con los parámetros que requiera en la sección de configuración del módulo en cuestión y habilitar el uso de Docker en las opciones de configuración de Nextflow.

Para descargar el docker que nos interesa
```
docker pull laugoro/resistance-workshop-inmegen:public
```

Para agregar la información en la sección de configuración del módulo, se debe agregar la siguiente línea en el archivo modules.nf en el proceso que requiere del docker
```
container 'laugoro/resistance-workshop-inmegen:public'
containerOptions "-v ${params.refdir}:/ref"
```

Para habilitar el uso de Docker en las opciones de configuración de Nextflow, se deben agregar las siguientes líneas en el archivo nextflow.config
```
docker {
    enabled = true
    temp = 'auto'
    fixOwnership = true
}
```


Visualiza el archivo modules.nf en la carpeta del ejercicio. Identifica las principales diferencias con respecto al ejercicio anterior:

1. Tiene dos procesos: fastp y bwa
2. El segundo proceso llama a un contenedor de Docker

```
// Align against genome

process bwa {
  cache 'lenient'
  container 'laugoro/resistance-workshop-inmegen:public'
  publishDir params.out, mode:'copy'

  input:
  tuple val(sample), path(fastp_data_R1), path(fastp_data_R2)

  output:
  tuple val(sample), path("alignment/${sample}_1.bam"),
                path("alignment/${sample}_1.bam.bai"), emit: bwa_out

  script:
  """

  mkdir -p alignment

  bwa mem  ${params.ref} \
        ${fastp_data_R1} \
        ${fastp_data_R2} \
        -t 12 |
  samtools sort -o alignment/${sample}_1.bam
  samtools index alignment/${sample}_1.bam alignment/${sample}_1.bam.bai
 """
}
```

¡Ahora si! Solo tienes que ejecutar el flujo de trabajo

```
nextflow run main.nf

```


## Flujo para la identificación de genes de resistencia

### Calidad y ensamblado


### Identificación de genes de resistencia


