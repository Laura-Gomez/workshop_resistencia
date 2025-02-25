# Taller ""Genómica aplicada para fortalecer la lucha contra las bacterias multirresistentes"

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
fastp \
        --in1 ${reads[0]} \
        --in2 ${reads[1]} \
        --out1 fastp_alone/${sample}_R1.fastq.gz \
        --out2 fastp_alone/${sample}_R2.fastq.gz \
        --html fastp_alone/${sample}.html \
        --json fastp_alone/${sample}.json
```

Ahora vamos a trasladar este comando a un flujo de trabajo de Nextflow, para esto necesitamos:

**1. Archivo main.nf**
Contiene los pasos del flujo de trabajo, el orden de los módulos que se van a ejecutar

**2. Archivo modules.nf**
Contiene las instrucciones específicas para cada módulo del flujo de trabajo

**3. Archivo nextflow.config**
Contiene los parámetros del flujo de trabajo

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




