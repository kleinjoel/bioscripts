#author: Joel Klein
#description: Used to indentify novel effector genes in translated RNA-seq sequences using Blast. Before running make sure seqeuence identifier in blast database does not contain |
# if it contains a | run ; eerst sed 's/cds.//' sed 's/|m.*|m./_/' velvetbremialetcontigs.fpa
#adds new table to annotation used the correct settings for signalp


 #!/usr/bin/perl
 #!/usr/local/bin/perl -w -I /Users/Joel/Perl5/

 #libraries
#add your custum BioPerl lib directory if neede
#use lib ("/home/../../perl5/Bio/SeqIO");


#updated version include tablemaking, annotation, and signalp annotation in table!

#no need to change
use warnings;
use strict;
#biotools used
#use Bio::Tools::Run::StandAloneBlastPlus;
use Bio::SearchIO;
use Bio::DB::Fasta;
use IO::String;
use Bio::SeqIO;
use Data::Dumper;

##start script
##declare variables
#related prot
#my 	$database = "all.fasta" ; #proteinseq
#input
#my $databaseold ="olddatabase"; #put here your sequences that you already identified
#my $databasePHI ="PHI/PHI"; #put here your sequences that you already identified
#make sure every file is present or comment all code related to a certain file to excluded it from the analysis
my $gff= "signalp5.gff3" ; # gff file with signalp info
my $wy= "pfs1_wy.tsv"; #table with WY sequences and locations
my $lwy= "LWY_all_proteins_pfs1.tsv"; #table with LWY sequences and locations

my $tmhmm="pfstmhmm2018.tsv"; #table with tmhmmhits
my $targetp="pfstargetp.tsv"; #table with targetphits
my $groupfile="Orthogroups.txt"; # file with orhomcl cluster data (space sep)

my $rxlrhmm="rxlr-hmm.tsv"; #list with genenames containing a hmm hit for rxlr
my $crnhmm="crn-hmm.tsv"; #list with genenames containing a hmm hit for crn

my $outdir = "outdir" ; #result files can be stored here
my $indir = "in" ; #input of sequences to be blasted
#translation hash!!
my $dict= "mrna_old-newnamesdict.tsv"; #table with WY sequences and locations
my $ANNOtbl= "Pfs1_annotation_information_table.tsv" ; # gff file with annotations info from blast



open( FILE1, "<", $dict ) or die "Can't find DICT table: $dict: $!\n";
##open( FILE2, "<", $blasttbl ) or die "Can't find DICT table: $blasttbl: $!\n";
my %hash;
my $key;
my $value;
while (<FILE1>) {
  chomp; # remove newline
  ($value, $key) = split '\t', $_, 2;
  $hash{$key} = $value;
}
#get annotation of table in dictornary

my $id;
my $sizecontig;
my $repeatregion;
my $combined;
my %hash2;
my $key2;
my $value2;
my $annotation;
my $shortname;

open( FILE2, "<", $ANNOtbl ) or die "Can't find Annotationtable table: $ANNOtbl: $!\n";

while (<FILE2>) {
   chomp;
      ## skipping the header lines
         my @ANNOtbl= split(/\t/);    ## splitting GFF line into array
         ## now use @gff as regular array
         #put each table entry in a different variable
         $id = $ANNOtbl[1];
         $sizecontig = $ANNOtbl[7];
         $repeatregion = $ANNOtbl[8];
         $shortname = $ANNOtbl[10];
         $annotation = $ANNOtbl[11];
         $combined = "$sizecontig\t$repeatregion\t$shortname\t$annotation";

         $hash2{$id} = $combined;
        # print "$name\n";
         #print "@blasttbl\n";
         #my @slice = @$ANNOtbl[7 .. 8];
         #print "@slice\n";
         ### specify the desired key

}#end of while for annotation table



##analises

##
print "BlastProgram is initilizing.. \nThe input direcotry is $indir \nCreated directory with found protein fasta files in $outdir" . "\nNow creating database from genome reads" . "\n";
#create a database from fasta file
#my $db = Bio::DB::Fasta->new($database);
#create table and add header
open (TABLEOUT, ">>$outdir/" . "table.tsv") or warn "could not create out file output file, $! \n";
print TABLEOUT "Query name\t Lenght\t Hashname\t Under 1kbp\t In repeat\t Shortname\t Annotation\t Group\t Signalp\t Cleavage site\t SPconf\t TargetP\t TPsignal\t WY-motif\t WY-number\t TMHMM\t Nhelix\t RXLR-motif\t RXLR-pos\t EER-motif\t EER-pos\t RXLRhmm\t CRN-motif\t CRN-pos\t HVL-motif\t HVL-pos\t CRNhmm\t KDEl\n";
opendir(DIR, $indir) or die "Cannot open directory to put the fasta sequences in $!" ;

##loop throught input folder and files to get sequences one by one
while (my $queryseq = readdir(DIR)) {
  # We only want files
   next unless (-f "$indir/$queryseq");
  # Use a regular expression to find files ending in .txt
   next unless ("$indir/$queryseq" =~ m/\.fasta$/);

  print "Input sequence is:\t$queryseq\n";

  my $queryseq2 = $queryseq ;
#    $queryseq2 =~ s/.fasta//;
  #go from the first aa seq to the next in a multiple fasta file query sequence
  my $seqio_obj = Bio::SeqIO->new(-file => "$indir/$queryseq", -format => "fasta" );
  while(my $seq_obj = $seqio_obj->next_seq){
    my $seqdesc= $seq_obj->desc ;
    my $seqid= $seq_obj->id ;
    my $seqstr1= $seq_obj->seq();
    my $seqlenght= length($seqstr1);
    ## do your  analisis
   print "Analizing:  ", "$seqid", "\t", $seqdesc, "\t";
   #my $score = "0";
   #my $crnscore = "0";

   print TABLEOUT "$seqid\t";
   print TABLEOUT "$seqlenght\t";
#get the hash name for the seqid
   if (exists($hash{$seqid}))
   {
     # if the key is found in the hash come here
     #print "Found Fred\n";
     print "$hash{$seqid}"."\n";
     print TABLEOUT "$hash{$seqid}\t";
     #print "@slice\n";
     #print "$keyid\n";
   }
   else
   {
     # come here if the key is not found in the hash
     print "Could not find Fred $seqid\n";
     print TABLEOUT "Could not find hash\t";
   }


   #get the annotation
   if (exists($hash2{$seqid}))
   {
     # if the key is found in the hash come here
     #print "Found Fred\n";
     print "$hash2{$seqid}"."\n";
     print TABLEOUT "$hash2{$seqid}\t";
     #print "@slice\n";
     #print "$keyid\n";
   }
   else
   {
     # come here if the key is not found in the hash
     print "Could not find Annotation $seqid\n";
     print TABLEOUT "Could not find Annotation\t";
   }


        #begin determining GROUPS
         open( GROUPS, "<", $groupfile ) or die "Can't find groupfile table: $groupfile: $!\n";

         my $groupname = "noGroup" ;

	 while (<GROUPS>) {
            chomp;
               ## skipping the header lines
		my @tmp = split(" ");


	foreach $a (@tmp){
		#print "$a ";
		my $gene_name = substr($a,0);
		#print "$seqid \t";
		#print "$gene_name \n";
		#print $gene_name ;
		if ($seqid eq $a ){
		#print $tmp[0] ." \n";
		print "OrthoMCL group found!\t";
		$tmp[0] =~ s/://g;
		$groupname="$tmp[0]";
	       } else {}

	} #end of foreach
      }
	print TABLEOUT "$groupname\t";
          # end of the GROUPS


      ##determine whether aa-seq has a signalpeptide using signalp
      my $signal = "no signal" ;
      my $cutoff = "0";
      my $cutoff2 = "0" ;
      my $cutoff3 = "0" ;
      my $conf = "0" ;

      #get the singalp annotation and put in tabel
      open( GFF, "<", $gff ) or die "Can't find gff singal p output: $gff: $!\n";
      while (<GFF>) {
         if (/^(\@)/) {
             chomp;
             }
         else {    ## skipping the header lines
               #s/##gff-version 2//;
               s/\n//;
               s/\r//;    ## removing new line
               my @gff = split(/\t+/);    ## splitting GFF line into array
               ## now use @gff as regular array
               #put each table entry in a different variable
               my $name     = $gff[0];
               my @namesplit = split(/_+/);
               my $nameshort = $namesplit[0];
               my $flag     = $gff[1];
               my $rname    = $gff[2];
               $conf    = $gff[5];
               $cutoff   = $gff[4];
               #print "$name"."\n";
               if ( $name eq $seqid ) {
                $signal = "Signal found" ;
                $cutoff2 = "$cutoff" ;    #cutoff in AA
                $cutoff3 = $cutoff2 * 3;  #cutoff in NT
                #$score = $score +5;
                #$crnscore = $crnscore +5;
               }
                #  else {
                ##print "NO signal found!!\n\n";
                # $cutoff2 = "0";
                #$cutoff3 = "0" ;
                #}

                 }
             }
            print TABLEOUT "$signal\t";
            print TABLEOUT "$cutoff2\t";
            print TABLEOUT "$conf\t";
            print "$signal \t";
            print "Cleavage site at $cutoff2 \n";
             # end of the determination signalp

         #begin determining TARGETP
         open( TARGETP, "<", $targetp ) or die "Can't find TARGETP table: $targetp: $!\n";
         my $nametarget;
         my $targetfound = "no TARGETP" ;
	 my $targetfound2 = "no TARGETP" ;
         my $signaltarget = "0";
	 my $signaltarget2 = "0";

	 while (<TARGETP>) {
            chomp;
               ## skipping the header lines
                  my @targetp = split(/\t/);    ## splitting GFF line into array
                  ## now use @gff as regular array
                  #put each table entry in a different variable
                  $nametarget =$targetp[0];
                  $targetfound =$targetp[6];
                  $signaltarget =$targetp[7];
                  #print "$namewy"."\t";
                  #print "$numberwy"."\n\n";
                  if ( $nametarget eq $seqid ) {
                   $targetfound2 = $targetfound;
                   $signaltarget2 = $signaltarget;

                   #print "$wyfound, $namewy, $numberwy2 \n";
                  }
                }
            print "targetp found: $nametarget $targetfound2\t";
            print "$signaltarget2\n";
            print TABLEOUT "$targetfound2\t";
            print TABLEOUT "$signaltarget2\t";

          # end of the determination TargetP

         #begin determining WY
         open( WY, "<", $wy ) or die "Can't find WY table: $wy: $!\n";
         my $namewy;
         my $wyfound = "no WY" ;
         my $numberwy = "0";
         my $numberwy2 = "0";
	 while (<WY>) {
            chomp;
               ## skipping the header lines
                  my @wy = split(/\t/);    ## splitting GFF line into array
                  ## now use @gff as regular array
                  #put each table entry in a different variable
                  $numberwy = $wy[0];
                  $namewy = $wy[1];
                  #print "$namewy"."\t";
                  #print "$numberwy"."\n\n";
                  if ( $namewy eq $seqid ) {
                   $wyfound = "WYfound" ;
                   $numberwy2 = "$numberwy";
                   #$score = $score +2;
                   print "$wyfound, $namewy, $numberwy2 \n";
                  }
                }
            print "$wyfound"."\t";
            print "$numberwy2"."\n";
            print TABLEOUT "$wyfound\t";
            print TABLEOUT "$numberwy2\t";

          # end of the determination WY
          #begin determining LWY
          open(LWY, "<", $lwy ) or die "Can't find LWY table: $lwy: $!\n";
          my $namelwy;
          my $lwyfound = "no LWY" ;
          my $numberlwy = "0";
          my $numberlwy2 = "0";
 	 while (<LWY>) {
             chomp;

                ## skipping the header lines
                   my @lwy = split(/\t/);    ## splitting GFF line into array
                   ## now use @gff as regular array
                   #put each table entry in a different variable
                   $numberlwy = $lwy[4];
                   $namelwy = $lwy[5];
						 #print "looking for LWY: $namelwy";
                   if ( $namelwy eq $seqid ) {
                    $lwyfound = "LWYfound" ;
                    $numberlwy2 = "$numberlwy";
                    #$score = $score +2;
                    print "$lwyfound, $namelwy, $numberlwy2 \n";
                   }
                 }
             print "$lwyfound"."\t";
             print "$numberlwy2"."\n";

             print TABLEOUT "$lwyfound\t";
             print TABLEOUT "$numberlwy2\t";

         # begin determination TMHMM
         open( TMHMM, "<", $tmhmm ) or die "Can't find WY table: $tmhmm: $!\n";
         my $nametmhmm;
         my $tmhmmfound = "no TMHMM" ;
         my $numberhelix = "0";
         my $numberhelix2 ="";
	 while (<TMHMM>) {
            chomp;
               ## skipping the header lines
                  my @tmhmm = split(/\t+/);    ## splitting GFF line into array
                  ## now use @gff as regular array
                  #put each table entry in a different variable
                  $nametmhmm = $tmhmm[0];
                  $numberhelix = $tmhmm[4];
                  #print "$namewy"."\t";
                  #print "$numberwy"."\n\n";
                  if ( $nametmhmm eq $seqid && $numberhelix>0 ) {
                   $tmhmmfound = "TMHMM found" ;
                   print "TMHMM found! N-helix: $numberhelix\n";
                   $numberhelix2 = $numberhelix;
                   #print "$wyfound, $namewy, $numberwy2 \n";
                  }
                }
            print "$tmhmmfound"."\t";
            print "$numberhelix2"."\n";
            print TABLEOUT "$tmhmmfound\t";
            print TABLEOUT "$numberhelix2\t";
          # end of the determination TMHMM



#done proccesing tables... now going on with the regex for rxlr,crn and kdel
      #	 #find RXLR regexs
         my $lenght=100;
         my $substr2=(substr($seqstr1,$cutoff2,$lenght));
         my $pos1 = "0";
         my $pos3 ="0";
#for RxLR

         if ($substr2 =~ m/..LR/) {
           #change the offset by changing the string with substring!
            print "match RXLR/GXLR\t";
            $pos1 = $-[0];
            print "$pos1\t";
            my $subrxlr = substr($substr2 ,$pos1 ,4);
            print "IN first: $subrxlr \n";
            $pos3=($pos1 + $cutoff2 +1);
            print TABLEOUT "$subrxlr\t" ;
            print TABLEOUT "$pos3 \t" ;
            #$score = $score +3;
          } elsif($substr2 =~ m/R.L./) {
            print "match RXLR/GXLR\t";
            $pos1 = $-[0];
            print "$pos1\t";
            my $subrxlr = substr($substr2 ,$pos1 ,4);
            print "IN SECOND: $subrxlr \n";
            $pos3=($pos1 + $cutoff2 +1);
            print TABLEOUT "$subrxlr\t" ;
            print TABLEOUT "$pos3 \t" ;
           # $score = $score +3;
         }else {
            print "no match RXLR/GXLR\t";
            print TABLEOUT "\t" ;
            print TABLEOUT "\t" ;
          }

         my $seqstr3=(substr($substr2,$pos1));
#for EER

            if ($seqstr3 =~ m/[ED][EDNQ][RK]/) {
           #change the offset by changing the string with substring!
            print "match EER like\t";
            my $pos5 = $-[0];
            my $pos6 =($pos5 + $cutoff2 + $pos1 +1);
            print "$pos6\t";
            my $subeer = substr($seqstr3 ,$pos5 ,3);
            print "$subeer \n";
            print TABLEOUT "$subeer\t" ;
            print TABLEOUT "$pos6 \t" ;
           # $score = $score +1;
          } else {
            print TABLEOUT "\t" ;
            print TABLEOUT "\t" ;
            print "no match EEX\n";
          }
            # begin determination RXLR HMMR
         open( RXLR, "<", $rxlrhmm ) or die "Can't find RXLRmmr table: $rxlrhmm: $!\n";
         my $namerxlr;
         my $rxlrfound = "no RXLRhmm" ;

	 while (<RXLR>) {
            chomp;
               ## skipping the header lines
                  my @rxlrhmm = split(/\t+/);    ## splitting GFF line into array
                  ## now use @gff as regular array
                  #put each table entry in a different variable
                  $namerxlr = $rxlrhmm[0];
                  if ( $namerxlr eq $seqid) {
                   $rxlrfound = "RXLRhmm found" ;
                  }
                }
            print "$rxlrfound"."\t";
            print TABLEOUT "$rxlrfound\t";


#for CRN
         my $pos4="0";
         if ($substr2 =~ m/.F.AK/) {
           #change the offset by changing the string with substring!
            print "match CRN-like\t";
            my $pos1 = $-[0];
            print "$pos1\t";
            my $subcrn = substr($substr2 ,$pos1 ,5);
            print "$subcrn \n";
           # $crnscore=$crnscore+5;
            $pos4=($pos1 + $cutoff2 +1);
            print TABLEOUT "$subcrn\t" ;
            print TABLEOUT "$pos4 \t" ;
          } else {
            print TABLEOUT "\t" ;
            print TABLEOUT "\t" ;
            print "no match LFLAK\n";
          }
         my $pos5="0";
         if ($substr2 =~ m/HVL/) {
           #change the offset by changing the string with substring!
            print "match CRN-like\t";
            my $pos1 = $-[0];
            print "$pos1\t";
            my $subdwl = substr($substr2 ,$pos1 ,3);
            print "$subdwl \n";
            $pos5=($pos1 + $cutoff2 +1);
            print TABLEOUT "$subdwl\t";
            print TABLEOUT "$pos5 \t";
           # $crnscore=$crnscore+2;
          } else {
            print "no matchHVL\n";
            print TABLEOUT "\t" ;
            print TABLEOUT "\t" ;
          }

          # begin determination RXLR HMMR
         open( CRN, "<", $crnhmm ) or die "Can't find CRNmmr table: $crnhmm: $!\n";
         my $namecrn;
         my $crnfound = "no CRNhmm" ;
	 while (<CRN>) {
            chomp;
               ## skipping the header lines
                  my @crnhmm = split(/\t+/);    ## splitting GFF line into array
                  ## now use @gff as regular array
                  #put each table entry in a different variable
                  $namecrn = $crnhmm[0];
                  if ( $namecrn eq $seqid) {
                   $crnfound = "CRnhmm found" ;
                  }
                }
            print "$crnfound"."\t";
            print TABLEOUT "$crnfound\t";

        # print TABLEOUT "$score\t" ;
        # print TABLEOUT "$crnscore\t" ;
        # $crnscore ="0";
        # $score = "0";

         my $pos6="0";



          # end of the determination RXLR
         #get KDEL motifs
         if ($seqstr1 =~ m/[KH]DEL$/) {
           #change the offset by changing the string with substring!
            print "match KDEL-motif\t";
            my $pos1 = $-[0];
            $pos6=($pos1 + $cutoff2 +1);
            print "$pos1\t";
            my $subkdel = substr($seqstr1 ,$pos1 ,4);
            print "$subkdel \n";
            print TABLEOUT "KDEL-motif: $pos1 \t" ;
          } else {
            print "no matchKDEL\n";
            print TABLEOUT "\t" ;
          }
#finished with Regex and HMM motifs for RxLR CRN and KDEL

#   #Figure out if the protein is already annotated
#    # process each seq
#      print "Figure out if the protein is already annotated:  ", "$seqid", "\t", $seqdesc, "\n";
#
#
#      #perform a blast search
#            my $percidenty2="1";
#	    my $percidenty3="0";
#            my $percidcheck="0";
#            my $iddesc="";
#            my $id="";
#            my $cloned="no";
#            my $hitlenght="0";
#            my $frac= "0";
#    my $result_old = $fac_old->blastp( -query => $seq_obj,
#				  -outfile => 'report_old.bls',
# 				  -method_args => [
# 				  '-evalue' => 1e-5,
#				  '-num_threads' => 4,
#                                  '-num_alignments' => 5]);
#
#    #print only significant hits
#    my  $report_obj2 = new Bio::SearchIO(-format => 'blast',
#				      -file   => 'report_old.bls');
#
#       #output parsing
#    while( $result_old = $report_obj2->next_result ) {
#	while(my $hit2 = $result_old->next_hit ) {
#	  while(my $hsp2 = $hit2->next_hsp ) {
##filtering
#	    if ($hsp2->evalue < 1e-5 && $percidcheck<$percidenty2) {
#	      #create output variables
#            $id = $hit2->name ;
#            $iddesc = $hit2->description;
#            $percidenty2 = $hsp2->percent_identity ;
#	    $percidenty3 = int($percidenty2) ;
#            $percidcheck = $percidenty2;
#            $hitlenght= $hit2->length ;
#            print  "\nold one found $percidenty3!!!\n";
#            }
#          }
#        }
#    }

###print already found in table!!####
#if ($percidenty2 > 80) {
#   $cloned = "Already cloned";
#            print TABLEOUT "$cloned\t";
#            print TABLEOUT "$id\t";
#            print TABLEOUT "$percidenty3\t";
#            $frac = $hitlenght/$seqlenght;
#            print TABLEOUT "$frac\t";
#   print  "\nAlready cloned!!\n";
#} else {
#   print TABLEOUT "not cloned\t";
#   print TABLEOUT "\t";
#   print TABLEOUT "\t";
#   print  "Unique and new!!!\n";
#}
#$cloned="";
#$percidenty3="";

#end of while loop!
#end of already cloned investigation!!!!

          ##END already cloned blast
          #start already PHI blast #http://www.phi-base.org/
#     print "Figure out if the protein is present in de PHIdatabse:"."$seqid"."\n\n";
#
#      #perform a blast search
#    my $result_PHI = $fac_PHI->blastp( -query => $seq_obj,
#				  -outfile => 'report_phi.bls',
# 				  -method_args => [
# 				  '-evalue' => 1e-15,
#				  '-num_threads' => 4,
#                                  '-num_alignments' => 3]);
#
#    #print only significant hits
#    my  $report_objPHI = new Bio::SearchIO(-format => 'blast',
#				      -file   => 'report_phi.bls');
#             my $idPHI;
#             my $percidenty;
#	     my $eval;
##
##       #output parsing
#    while( $result_PHI =  $report_objPHI->next_result ) {
#	while(my $hit3 = $result_PHI->next_hit ) {
#	  while(my $hsp3 = $hit3->next_hsp ) {
#
#
###filtering hits
#	    if ($hsp3->evalue < 1e-10) {
#	      #create output variables
#             $idPHI = $hit3->name ;
#             $percidenty = $hsp3->percent_identity ;
#	     $eval = $hsp3->evalue ;
#             print TABLEOUT "$idPHI\t";
#             print TABLEOUT "$eval\t";
#            }
#          }
#        }
#    }

   ##end of analis
   print TABLEOUT "\n";
  }#end of sequence process loop
  }#end of sequence input loop
#$fac_old->cleanup;
#$fac_PHI->cleanup;
close (TABLEOUT);
close (GFF);
close FILE1;
close FILE2;
close FILE3;
close (WY);
close (RXLR);
close (CRN);
closedir(DIR);
print "fineshed closing program...\n";
sleep(3);
exit;
