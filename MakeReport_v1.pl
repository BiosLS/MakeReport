#!/usr/bin/perl -w
use strict;
##############################################
#Desscription:
Generates a frequency table of affected genes on ChAS 
generated segment table data. Data is derived from 
affymetrix's Cytoscan HD microarray
##############################################
#Nombre del script= MakeGlobalReport_v1.pl
#Date: July 20 2020
#Version= 1.0
##############################################
#Usage
##############################################
my $num_args = $#ARGV + 1;
if ($num_args != 1) {
    print "\nUsage:\n\n";
    print "perl MakeGlobalReport_v1.pl Test\n\n";
    print "perl MakeReport.pl - name of script\n";
    print "Test - Folder with segment tables in txt format generated with ChAS\n";
    exit;
}
my $carpeta = $ARGV[0];

my $Separador= ",";

##############################################
#Database
##############################################
my $DechipherFileBed = "HI_Predictions_Version3.bed";
my $sourceDEchipher = "Databases/$DechipherFileBed";
open(FILE1, "<$sourceDEchipher" ) or die( "Could not open $sourceDEchipher: $!\n" );
#open(FILE1, $DechipherFileBed);
my @Archivobed = <FILE1>;
chomp @Archivobed;
my $encabezado2 = shift @Archivobed;

my $gdcfile = "GastricAdenocarcinoma_portal.gdc.cancer.tsv";
my $sourcegdcfile = "Databases/$gdcfile";
open(FILE2, "<$sourcegdcfile" ) or die( "Could not open $sourcegdcfile: $!\n" );
#open(FILE2, $gdcfile);
my @DataBasegdcfile = <FILE2>;
chomp @DataBasegdcfile;

my $gdcDLBCfile = "emt.hsa.info_seq.txt";
my $sourcegdcDLBCfile = "Databases/$gdcDLBCfile";
open(FILE3, "<$sourcegdcDLBCfile" ) or die( "Could not open $sourcegdcDLBCfile: $!\n" );
my @DataBasegdcDLBCfile = <FILE3>;
chomp @DataBasegdcDLBCfile;

#my $DDG2Pfile = "DDG2P_9_8_2017.csv";
#my $sourceDDG2Pfile = "Databases/$DDG2Pfile";
#open(FILE4, "<$sourcegdcDLBCfile" ) or die( "Could not open $sourceDDG2Pfile: $!\n" );
#my @DataBaseDDG2Pfile = <FILE4>;
#chomp @DataBaseDDG2Pfile;

my $NumeroArchivos=0;

my $Encabezado1 = "Gene".$Separador."Frequency".$Separador."%(n=".$NumeroArchivos.")".$Separador."Chromosome".$Separador."Cytoband".$Separador."Type".$Separador."State".$Separador."OMIM phenotype".$Separador."Haploinsufficiency Predictions".$Separador."Gene Name".$Separador."Affected Gastric Cases in Cohort GDC".$Separador."Affected Gastric Cases Across GDC".$Separador."MesenchimalGene";
##############################################
##############################################
#CV analysis
##############################################
####Stage 1: Report of CNV###########
print "Stage 1: Report of CNV\n";
my $CNVfile = "$carpeta/cnv.txt";
open(FILE, "<$CNVfile") or die( "No se encuentra el $CNVfile archivo: $!\n" );
my @ArchivoSegmentos = <FILE>;

my $contadorSegNC = scalar @ArchivoSegmentos;

my $Gain_CNS_Total;
my $Loss_CNS_Total;
my $ArchivosEstudiados;

my $OuputFileName = "$carpeta/Report_TableFrequency_CNV.csv";

if($contadorSegNC < 1){
    print "File cnv.txt with cero events\n";
    $NumeroArchivos = 0;
    $Gain_CNS_Total=0;
    $Loss_CNS_Total=0;
    #open (OUTPUTFILE, ">$OuputFileName");
    open(FILE, ">$OuputFileName");
    print FILE "File mosaic.txt with cero events\n";
}

if($contadorSegNC >= 1){
my $encabezado = shift @ArchivoSegmentos;
$contadorSegNC=$contadorSegNC-1;
my @Headear = BuscaHeaders ($encabezado);
#Encabezdos de interes
#$File=$Headear[0];
#$Genes=$Headear[1];
#$CN_State=$Headear[2];
#$Type=$Headear[3];
#$Chromosome=$Headear[4];
#$Cytoband=$Headear[5];
#$OMIMphenotipe=$Headear[6];

my @GenesUnicos = subrutinaGenesUnicos ($Headear[0], $Headear[1], $Headear[2], $Headear[3], $Headear[4], $Headear[5], $Headear[6], @ArchivoSegmentos);
my $outputfileGU = "genesunics.txt";
open (OUTPUT1, ">$outputfileGU");
foreach (@GenesUnicos) {
  print OUTPUT1 "$_\n";
}

my @info = subrutinaNoFilesNoeventos ($contadorSegNC, $Headear[0], $Headear[3], \@ArchivoSegmentos);
$NumeroArchivos = $info[0];
$Gain_CNS_Total = $info[1];
$Loss_CNS_Total = $info[2];
$ArchivosEstudiados = $info[6];
#print "numero de archivos $info[0]\n";

my $Encabezado1 = "Gene".$Separador."Frequency".$Separador."%(n=".$NumeroArchivos.")".$Separador."Chromosome".$Separador."Cytoband".$Separador."Type".$Separador."State".$Separador."OMIM phenotype".$Separador."Haploinsufficiency Predictions".$Separador."Gene Name".$Separador."Affected gastric Cases in Cohort GDC".$Separador."Affected gastric Cases Across GDC".$Separador."Gene Mesenchimal".$Separador."Info";
my @TablaFrecuencias = subrutinaTablaFrecuencias ($Separador, $info[0], $OuputFileName, $Encabezado1, \@GenesUnicos);
close FILE;
}

######################################################################################
print "Paso 2: Generando Reporte de frecuencias de Mosaicismo\n";
######################################################################################

my $Mosaicfile = "$carpeta/mosaic.txt";

open(FILE, "<$Mosaicfile") or die( "No se encuentra el $Mosaicfile archivo: $!\n" );
@ArchivoSegmentos = <FILE>;
$contadorSegNC = scalar @ArchivoSegmentos;
#print "$contadorSegNC\n";

$OuputFileName = "$carpeta/Report_TablaFrequency_Mosaic.csv";

my $Gain_Mosaic_Total;
my $Loss_Mosaic_Total;

if($contadorSegNC < 1){
    print "File mosaic.txt cero events\n";
    $NumeroArchivos = 0;
    $Gain_Mosaic_Total = 0;
    $Loss_Mosaic_Total = 0;
    #open (OUTPUTFILE, ">$OuputFileName");
    open(FILE, ">$OuputFileName");
    print FILE "File mosaic.txt cero events\n";
}

#if($contadorSegNC >= 1){
#my $encabezado = shift @ArchivoSegmentos;
#$contadorSegNC = $contadorSegNC - 1;
#my @Headear = BuscaHeaders ($encabezado);
#my @GenesUnicos = subrutinaGenesUnicos ($Headear[0], $Headear[1], $Headear[2], $Headear[3], $Headear[4], $Headear[5], $Headear[6], @ArchivoSegmentos);

#my @info = subrutinaNoFilesNoeventos ($contadorSegNC, $Headear[0], $Headear[3], \@ArchivoSegmentos);
#$NumeroArchivos = $info[0];
#$Gain_Mosaic_Total = $info[3];
#$Loss_Mosaic_Total = $info[4];
#my $Encabezado1 = "Gene".$Separador."Frequency".$Separador."%(n=".$NumeroArchivos.")".$Separador."Chromosome".$Separador."Cytoband".$Separador."Type".$Separador."State".$Separador."OMIM phenotype".$Separador."Haploinsufficiency Predictions".$Separador."Gene Name".$Separador."Affected AML Cases in Cohort GDC".$Separador."Affected AML Cases Across GDC".$Separador."Gene Mesenchimal".$Separador."info";
#my @TablaFrecuencias = subrutinaTablaFrecuencias($Separador, $info[0], $OuputFileName, $Encabezado1, \@GenesUnicos);
#close FILE;
#}
##############################################
#LOH analysis
##############################################
print "Stage 3: Generating Report of frequencies of Loss of Heterogocity\n";

#open(FILE, "TablesegmentsLOH.txt");
my $LOHfile = "$carpeta/loh.txt";
open(FILE, "<$LOHfile") or die( "No se encuentra el $LOHfile archivo: $!\n" );
@ArchivoSegmentos = <FILE>;
$contadorSegNC = scalar @ArchivoSegmentos;
my $LOH_Total;

#if($contadorSegNC < 1){
#    print "File loh.txt cero events\n";
#    $NumeroArchivos = 0;
#    $LOH_Total = 0;
    #open (OUTPUTFILE, ">$OuputFileName");
#    open(FILE, ">$OuputFileName");
#    print FILE "El archivo loh.txt no tiene eventos\n";
#}

#if($contadorSegNC >= 1){

#my $encabezado = shift @ArchivoSegmentos;
#$contadorSegNC = $contadorSegNC-1;

#my @Headear = BuscaHeaders ($encabezado);

#$OuputFileName = "$carpeta/Reporte_TablaFrecuencias_LOH.csv";
#my @GenesUnicos = subrutinaGenesUnicos ($Headear[0], $Headear[1], $Headear[2], $Headear[3], $Headear[4], $Headear[5], $Headear[6], @ArchivoSegmentos);

#my @info = subrutinaNoFilesNoeventos ($contadorSegNC, $Headear[0], $Headear[3], \@ArchivoSegmentos);
#$NumeroArchivos = $info[0];
#$LOH_Total = $info[5];
#my $Homozygous_deletions=0;
#my $Heterozygous_deletions=0;
#my $CNN_LOH=0;
#my $Encabezado1 = "Gene".$Separador."Frequency".$Separador."%(n=".$NumeroArchivos.")".$Separador."Chromosome".$Separador."Cytoband".$Separador."Type".$Separador."State".$Separador."OMIM phenotype".$Separador."Haploinsufficiency Predictions".$Separador."Gene Name".$Separador."Affected gastric Cases in Cohort GDC".$Separador."Affected gastric Cases Across GDC".$Separador."gene mesenchimal".$Separador."info";
#my @TablaFrecuencias = subrutinaTablaFrecuencias($Separador, $info[0], $OuputFileName, $Encabezado1, \@GenesUnicos);
#close FILE;
#}

##############################################
#Report global
##############################################
print "Stage 4: Make Global Report\n";

my $outputfile = "$carpeta/Reporte_global.csv";
open (OUTPUT1, ">$outputfile");

#Encabezado 1
$Encabezado1 = "All ".$NumeroArchivos." files";
#Encabezado 1
my $Encabezado2 = "Abnormalites".$Separador."N".$Separador."Per Sample";
my $CNS_Total= $Gain_CNS_Total + $Loss_CNS_Total + $Gain_Mosaic_Total + $Loss_Mosaic_Total + $LOH_Total;
my $PerSample = 1;

#Lineas Reporte

my $PerSampleCNGains = $Gain_CNS_Total/$CNS_Total;
$PerSampleCNGains = sprintf("%.2f", $PerSampleCNGains);
my $lineCNVGain = "Totals CNV Gains".$Separador.$Gain_CNS_Total.$Separador.$PerSampleCNGains;

my $PerSampleCNloss = $Loss_CNS_Total/$CNS_Total;
$PerSampleCNloss = sprintf("%.2f", $PerSampleCNloss);
my $lineCNVloss = "Totals CNV Losses".$Separador.$Loss_CNS_Total.$Separador.$PerSampleCNloss;

my $PerSampleMosaicGain = $Gain_Mosaic_Total/$CNS_Total;
$PerSampleMosaicGain = sprintf("%.2f",$PerSampleMosaicGain);
my $lineMosaicGain = "Totals Mosaic Gains".$Separador.$Gain_Mosaic_Total.$Separador.$PerSampleMosaicGain;

my $PerSampleMosaicLoss = $Loss_Mosaic_Total/$CNS_Total;
$PerSampleMosaicLoss = sprintf("%.2f",$PerSampleMosaicLoss);
my $lineMosaicLoss = "Totals Mosaic Losses".$Separador.$Loss_Mosaic_Total.$Separador.$PerSampleMosaicLoss;

my $PerSampleLOH = $LOH_Total/$CNS_Total;
$PerSampleLOH = sprintf("%.2f",$PerSampleLOH);
my $lineCNNLOH = "LOH".$Separador.$LOH_Total.$Separador.$PerSampleLOH;
my $lineTotal = "Total".$Separador.$CNS_Total.$Separador.$PerSample;
my $linesFilesEstudiado = "Files".$Separador.$ArchivosEstudiados;

print OUTPUT1 "$Encabezado1\n$Encabezado2\n$lineCNVGain\n$lineCNVloss\n$lineMosaicGain\n$lineMosaicLoss\n$lineCNNLOH\n$lineTotal\n\n$linesFilesEstudiado";

close OUTPUT1;

exit;

######################################################################
######################################################################
######Subrutins
######################################################################
######################################################################

################################################
#Localiza la posicion de las columnas de interes
################################################
sub BuscaHeaders{
my $encabezado = shift @_;
   
  my @segmentosHeader= split ("\t", $encabezado);
  chomp  @segmentosHeader;
  my $contador1 = scalar @segmentosHeader;
    #columnas de interes de la tabla de segmentos
    my $File=0;
    my $Genes=0;
    my $CN_State=0;
    my $Type=0;
    my $Chromosome=0;
    my $Cytoband=0;
    my $OMIMphenotipe=0;
    
    for (my $i=1; $i < $contador1; $i++){
        my $columnas = $segmentosHeader[$i];
        chomp $columnas;
        if($columnas eq 'File'){$File = $i;}
        if($columnas eq 'CN State'){$CN_State = $i;}
        if($columnas =~ 'Type'){$Type = $i;}
        if($columnas =~ "Chromosome"){$Chromosome = $i;}
        if($columnas =~ "OMIM" and $columnas =~ "Phenotype"){$OMIMphenotipe = $i;}
        #if($columnas eq 'OMIM'){$OMIMphenotipe = $i;}
        if(($columnas eq 'Genes') and ($columnas ne 'OMIM') and ($columnas ne 'Count') ){$Genes = $i;}
        if($columnas =~ "Cytoband"){$Cytoband = $i;} 
    }
    
    my @headers = ("$File","$Genes","$CN_State","$Type","$Chromosome","$Cytoband","$OMIMphenotipe");
    return @headers;
}#subrutina

################################################
#Genera La lista de Genes Unicos
################################################
sub subrutinaGenesUnicos {
    my( $File, $Genes, $CN_State, $Type, $Chromosome, $Cytoband, $OMIMphenotipe, @ArchivoSegmentos) = @_;
#$File=$Headear[0];
#$Genes=$Headear[1];
#$CN_State=$Headear[2];
#$Type=$Headear[3];
#$Chromosome=$Headear[4];
#$Cytoband=$Headear[5];
#$OMIMphenotipe=$Headear[6];
    my $joiner = ",";
    my @ListaGenesUnicos;
    my $contador2 = scalar @ArchivoSegmentos;
    
    for (my $ii=0; $ii < $contador2; $ii++){
        my $renglon = $ArchivoSegmentos[$ii];
        my @segmentosBody = split ("\t", $renglon);
        chomp  @segmentosBody;
        my $archivo = $segmentosBody[$File];
        $archivo =~ s/\s//g;
        my $genestodos = $segmentosBody[$Genes];
        my $CN_estado = $segmentosBody[$CN_State];
        $CN_estado =~ s/\s//g;
        my $tipodevento = $segmentosBody[$Type];
        $tipodevento =~ s/\s//g;
        my $Chromosoma = $segmentosBody[$Chromosome];
        $Chromosoma =~ s/\s//g;
        my $Cytobanda = $segmentosBody[$Cytoband];
        $Cytobanda=~ s/\s//g;
        my $OMIMphenotipe_info  = $segmentosBody[$OMIMphenotipe];
        $OMIMphenotipe_info  =~ s/\([^)]+\)//g;
        $OMIMphenotipe_info  =~ s/, \{/|/g;
        $OMIMphenotipe_info  =~ s/\} ,/|/g;
        $OMIMphenotipe_info  =~ s/\}//g;
        $OMIMphenotipe_info  =~ s/\] ,/|/g;
        $OMIMphenotipe_info  =~ s/, \[/|/g;
        $OMIMphenotipe_info  =~ s/\{/|/g;
        $OMIMphenotipe_info  =~ s/\}/|/g;
        $OMIMphenotipe_info  =~ s/\[/|/g;
        $OMIMphenotipe_info  =~ s/\]/|/g;
        $OMIMphenotipe_info  =~ s/,//g;  
        $OMIMphenotipe_info  =~ s/\?//g;    
        my @genessingle = split (",", $genestodos);
        chomp  @genessingle;
        my $contador3 = scalar @genessingle; 
        for (my $iii=0; $iii < $contador3; $iii++){
            $genessingle[$iii] =~ s/\s//g;
            #my $fullNIPname = $genessingle[$iii].$joiner.$Chromosoma.$joiner.$Cytobanda;
            my $fullNIPname = $genessingle[$iii].$joiner.$Chromosoma.$joiner.$Cytobanda.$joiner.$tipodevento.$joiner.$CN_estado.$joiner.$OMIMphenotipe_info.$joiner.$archivo;
            push (@ListaGenesUnicos, $fullNIPname);
        }
                  
    }
    print FILE2 "@ListaGenesUnicos\n";
    return @ListaGenesUnicos;
}

########################################
#subrutina para tomar datos: Numero de archivos(muestras), No eventos perdida, No eventos ganancia
########################################

sub subrutinaNoFilesNoeventos {
  my($contadorSegNC, $File, $Type) = @_;
  my @ArchivoSegmentos = @{$_[3]};
  my @info;

#conteo de typo de evento
  my @filestodos;
  my $Gain_CNS_Total = 0;
  my $Loss_CNS_Total = 0;
  my $GainMosaic_Total= 0;
  my $LossMosaic_Total= 0;
  my $LOH_Total = 0;
  
  for (my $ii=0; $ii < $contadorSegNC; $ii++){
    my $renglon = $ArchivoSegmentos[$ii];
    my @segmentosBody = split ("\t", $renglon);
    chomp  @segmentosBody;
    my $Type_info = $segmentosBody[$Type];
    $Type_info =~ s/\s//g;
    if ($Type_info  =~ "Loss"){$Gain_CNS_Total++;}
    if ($Type_info  =~ "Gain"){$Loss_CNS_Total++;}
    if ($Type_info  =~ "GainMosaic"){$GainMosaic_Total++;}
    if ($Type_info  =~ "LossMosaic"){$LossMosaic_Total++;}
    if ($Type_info  =~ "LOH"){$LOH_Total++;}
    push (@filestodos, $segmentosBody[$File]);
  }
  #eliminando archivos
  my %seen;
  @filestodos = grep { ! $seen{ $_ }++ } @filestodos;
  chomp (@filestodos);
  my $NumeroArchivos = scalar @filestodos;
  my $filejoined = join (';',@filestodos);
  #print "Numero de muestras analizadas: $NumeroArchivos\n";
  @info = ("$NumeroArchivos", "$Gain_CNS_Total", "$Loss_CNS_Total", "$GainMosaic_Total", "$LossMosaic_Total", "$LOH_Total", "$filejoined");
  return @info;
}#subrutina

################################################
#Tabla de Frecuencias
################################################

sub subrutinaTablaFrecuencias {

    my ($Separador, $NumeroArchivos, $OuputFileName, $Encabezado1) = @_;

    my @GenesUnicos = @{$_[4]};
    chomp @GenesUnicos;
    open (OUTPUTFILE, ">$OuputFileName");
    print OUTPUTFILE "$Encabezado1\n";
  
    my $sep = ",";
    my @archivos = ();
    my $tipos = "";
    my $states = "";
     my @final = ();
    foreach  my $element (@GenesUnicos) {
        #my $element= shift @GenesUnicos;
        @archivos = ();
        $tipos = "";
        $states = "";
        my @segmentos = split (",",  $element);
        chomp  @segmentos;
        my $baseName = $segmentos[0].$segmentos[1].$segmentos[2];
        $segmentos[6] =~ s/\s//g;
        my $OMi1 = $segmentos[6];
        push @archivos, $segmentos[6];
    
        foreach my $element2 (@GenesUnicos) {
            #my $element2= shift @GenesUnicos;;
            my @segmentos2 = split (",",  $element2);
            chomp  @segmentos2;
            my $baseName2 = $segmentos2[0].$segmentos2[1].$segmentos2[2];
            my $OMi2 = $segmentos[6];
        
            if ($baseName eq $baseName2 and $OMi1 eq $OMi2){
                $segmentos2[6] =~ s/\s//g;
                push @archivos, $segmentos2[6];
                $tipos = $tipos."|".$segmentos2[3];
                $states = $states."|".$segmentos2[4];
                next;
            }
            #if ($baseName ne $baseName2){
             #   push @GenesUnicos, $element2;
            #}
        
        }
        my %seen;
        @archivos = grep { ! $seen{ $_ }++ } @archivos;
        chomp (@archivos);
        my $contador = scalar @archivos;
    
        my $frecuencia = ($contador*100)/$NumeroArchivos;
        $frecuencia = sprintf("%.1f", $frecuencia);
        
        ##########################################################################################
        #Agregando Info de INformacion de la base de datos de DECIPHER de Haploincuficeincia
        ##########################################################################################
        
        my $HaploInfo = '';
        foreach my $fila2(@Archivobed) {
            my @Bedcontenido = split ("\t",  $fila2);
            chomp  @Bedcontenido;
            
            my @Probabilidades = split (/\|/,  $Bedcontenido[3]);
            chomp @Probabilidades;
    
            if ($segmentos[0] eq $Probabilidades[0]){
                $HaploInfo = $Bedcontenido[3];
                last;
            }
            if(  (\$fila2 == \$Archivobed[-1]) and ($segmentos[0] ne $Probabilidades[0]) ) {
                $HaploInfo = "NA";
                
            }
        }
        ###############################################################################
        #Anotando datos de Harmonized Cancer Datasets, Genomic Data Commons Data Portal
        #National Cancer Institude
        ##############################################################################
        my $DescripcionGene = '';
        my $Percent = '';
        my $Cases = '';
        
          foreach my $fila2(@DataBasegdcfile) {
            my @gdccontenido = split ("\t",  $fila2);
            chomp  @gdccontenido;
            
            if ($segmentos[0] eq $gdccontenido[0]){
                $gdccontenido[1] =~ s/,//g;
                $DescripcionGene = $gdccontenido[1];
                $Percent = $gdccontenido[4];
                $gdccontenido[5] =~ s/,//g;
                $Cases = $gdccontenido[5];
                last;
            }
            if(  (\$fila2 == \$DataBasegdcfile[-1]) and ($segmentos[0] ne $gdccontenido[0]) ) {
                $DescripcionGene = 'NA';
                $Percent = 'NA';
                $Cases = 'NA';
                
            }
        }
        #Lymphoid Neoplasm Diffuse Large B-cell Lymphoma
        my $PercentDLBC = '';
        my $CasesDLBC = '';
        
        foreach my $fila2(@DataBasegdcDLBCfile) {
            my @gdcDLBCcontenido = split ("\t",  $fila2);
            chomp  @gdcDLBCcontenido;
            $gdcDLBCcontenido[1] =~ s/\s//g;
            if ($segmentos[0] eq $gdcDLBCcontenido[1]){
                $gdcDLBCcontenido[1] =~ s/,//g;
                #$DescripcionGene = $gdcDLBCcontenido[1];
                $PercentDLBC = $gdcDLBCcontenido[1];
                $gdcDLBCcontenido[12] =~ s/,//g;
                $CasesDLBC = $gdcDLBCcontenido[12];
                last;
            }
            if(  (\$fila2 == \$DataBasegdcDLBCfile[-1]) and ($segmentos[0] ne $gdcDLBCcontenido[0]) ) {
                #$DescripcionGene = 'NA';
                $PercentDLBC = 'NA';
                $CasesDLBC = 'NA';
                
            }
        }
        
        my $diseasename = '';
        my $DDDcategory = '';
        my $organspecificity = '';
        
        ####Imprimiendo######
         
        my $lineaprint = $segmentos[0].$sep.$contador.$sep.$frecuencia.$sep.$segmentos[1].$sep.$segmentos[2].$sep.$tipos.$sep.$states.$sep.$segmentos[5].$sep.$HaploInfo.$sep.$DescripcionGene.$sep.$Percent.$sep.$Cases.$sep.$PercentDLBC.$sep.$CasesDLBC.$sep.$diseasename.$sep.$DDDcategory.$sep.$organspecificity;
        push( @final , $lineaprint);
        #print OUTPUTFILE "$lineaprint\n";
        @archivos = ();
    }
            
   my %seen;
  my @finalfinal = grep { ! $seen{ $_ }++ } @final;
  chomp (@finalfinal);
  
  foreach (@finalfinal) {
  print OUTPUTFILE "$_\n";
}
  
close OUTPUTFILE;

}#subrutina


