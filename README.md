# Antigenic-Cartography-WDL
:star: Parallel Antgenic Cartography with QC checks orchestrated by WDL :star:
## I. INTRODUCTION
This pipeline uses one R Version 4.3.1 [docker image](https://github.com/Cameron-Nguyen1/Docker_Images/tree/main/WDL_R_AntigenicCartography) to run an antigenic cartography analysis.</br>
</br>
An input titer table is specified, alongside cartography parameters, to derive a cartography figure,</br>
within group distances, and associated QC checks.</br>
</br>
The script is designed to run in parallel though the QC step is only bottlenecked by the dimension test.</br>

## II. INSTRUCTIONS:
Pass input titer table(s) to the script WDL script by using the "Input_AC.json" file.</br>
Use a WDL engine, like miniWDL, to run the workflow like so:</br>
`miniwdl run AntigenicCartography.wdl --input Input_AC.json`

## III. INPUT ARGUMENTS:
```
{
    "agQC.titertables": [LIST_OF_INPUT(S)],
    "agQC.prefixes": [LIST_OF_UNIQUE_IDENTIFIERS],
    "agQC.cpu": "1",
    "agQC.mem": "2G",
    "agQC.repsxdim": "100",
    "agQC.xy_lim": "-10,10,-10,10",
    "agQC.psizes": "5,2",
    "agQC.opacity": ".8,1",
    "agQC.agoverprint": "FALSE",
    "agQC.agsort": "TRUE"
}
```
Where:</br>
</br>
`titertables` : will be a list of inputs. </br>
`prefixes` : are unique identifiers respective to titer tables. </br>
`cpu` : number of CPUs to allocate to any individual task. </br>
`mem` : memory in GB to allocate to any individual task. </br>
`repsxdim` : number of test replicates per dimension in QC task, default should be 1000. </br>
`xy_lim` : XY limits for the cartography map in the following format: -X,X,-Y,Y </br>
`psizes` : Antigenic and Sera point sizes in the format: AntigenSize,SeraSize </br>
`opacity` : Controls opacity of Antigen and Sera poins in the format: AntigenOpacity,SeraOpacity </br>
`agoverprint` : Controls Antigen overprint, TRUE if Antigen should print over Sera, FALSE if not.</br>
`agsort` : Controls whether Antigens will be sorted by name. TRUE/FALSE. </br> 

## IV. EXAMPLE USAGE:
```
miniwdl run AntigenicCartography.wdl --input Input_AC.json
```

