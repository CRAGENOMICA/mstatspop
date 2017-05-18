#!/usr/local/bin/perl
use Getopt::Long;
use File::Find;
use strict;

my (%options,@colnames,@colnames_real,$ncol,@coldata,$nrow,@colnumbers);

#options
GetOptions(
     'in=s'      => \$options{'in'},
     'fc=s'      => \$options{'fc'},
);

#Usage
if($options{'in'} eq undef || $options{'fc'} eq undef) {
   print "\nUsage: \nperl collect_data_columns.pl -in (inputfile in mstatspop format)  -fc (file with the text of each header column each row) \n\n";
   exit;
}

#read the column names file
$ncol = 0;
open(INPUTC,"$options{'fc'}") or die "\nCan't open ".$options{'fc'}." file.\n";
while(<INPUTC>) {
   chomp($_);
   push @colnames,$_;
   $ncol = $ncol + 1;
}
close(INPUTC);


#read the mstatspop and keep the necessary columns
$nrow=0;
open(INPUT,"$options{'in'}") or die "\nCan't open ".$options{'in'}." file.\n";
while(<INPUT>) {
   my $line = $_;
   chomp($line);
   my @col_lines = split /\t/,$line;
   if($nrow==0) {
      for(my $col=0;$col<$ncol;$col++) {
         for(my $colv=0; $colv<=$#col_lines;$colv+=2) {
            if($col_lines[$colv] =~ /\Q$colnames[$col]\E/) {#eq $colnames[$col]) {#
               push @colnumbers, $colv+1;
               push @colnames_real, $col_lines[$colv];
            } 
         }  
      }
      $ncol = $#colnumbers + 1;
   }
   my @coldata_lines;
   foreach my $colv (@colnumbers) {
      push @coldata_lines, $col_lines[$colv]
   }
   push @coldata, @coldata_lines;
   $nrow = $nrow + 1;
}
close(INPUT);

#print selected data
for(my $colv=0;$colv<$ncol;$colv++){
   print $colnames_real[$colv]."\t";
}
print "\n";
for(my $row=0;$row<$nrow;$row++) {
   for(my $colv=0;$colv<$ncol;$colv++){
      print $coldata[$row*$ncol+$colv]."\t";
   }
   print "\n";
}

exit(0);
