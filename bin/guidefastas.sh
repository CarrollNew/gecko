#!/bin/bash 

BINDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

if [ $# -lt 7 ]; then
   echo " ==== ERROR ... you called this script inappropriately."
   echo ""
   echo "   usage:  $0 seqXName.fasta seqYName.fasta guideFile dimension LEN SIM WL"
   echo ""
   exit -1
fi






seqNameX=$(basename "$1")
extension="${seqNameX##*.}"
seqNameX="${seqNameX%.*}"


sed '/>/d' $1 > $seqNameX.temp
tr -d '\n' < $seqNameX.temp > $seqNameX.fix

seqNameY=$(basename "$2")
seqNameY="${seqNameY%.*}"


sed '/>/d' $2 > $seqNameY.temp
tr -d '\n' < $seqNameY.temp > $seqNameY.fix

rm $seqNameX.temp $seqNameY.temp


guided=$3
dimension=$4

LEN=$5
SIM=$6
WL=$7

lenX=$(wc -c $seqNameX.fix | awk '{print $1}')
lenY=$(wc -c $seqNameY.fix | awk '{print $1}')

echo "Length X: $lenX, Length Y: $lenY"

counterX=0
counterY=0

turn=0

mkdir tempfastas
mkdir all-results

echo "All by-Identity Ungapped Fragments (Hits based approach)
[Abr.98/Apr.2010/Dec.2011 -- <ortrelles@uma.es>
SeqX filename        : undef
SeqY filename        : undef
SeqX name            : undef
SeqY name            : undef
SeqX length          : $lenX
SeqY length          : $lenY
Min.fragment.length  : undef
Min.Identity         : undef
Tot Hits (seeds)     : undef
Tot Hits (seeds) used: undef
Total fragments      : undef
========================================================
Total CSB: 0
========================================================
Type,xStart,yStart,xEnd,yEnd,strand(f/r),block,length,score,ident,similarity,%ident,SeqX,SeqY" > all-results/master.csv

ratioX=$(($lenX / $dimension))
ratioY=$(($lenY / $dimension))
actualX=0
actualY=0
echo "Ratio X: $ratioX, Ratio Y: $ratioY"

for i in $( tail -n +2 $guided ); do
	#echo "$i"

	if [[ turn -eq 0 ]]; then
	
		actualX=`expr $i - 1` 
		actualX=$(($actualX * $ratioX))

	
		echo ">nothing" > tempfastas/X_${counterX}.fasta
		(tail -c +"$actualX" "$seqNameX.fix" | head -c "$ratioX") >> tempfastas/X_${counterX}.fasta
		
	
		counterX=`expr $counterX + 1`
	fi
	
	if [[ turn -eq 1 ]]; then
		
		actualY=`expr $i - 2` # used to be -2
		actualY=$(($actualY * $ratioY))
		fakeActualY=`expr $i + 2`
		fakeY=$(($fakeActualY * $ratioY))
		
		dasSUMX=`expr $actualX + $ratioX`
		dasSUMY=`expr $actualY + $ratioY`
		echo "coord X,Y: $actualX-$dasSUMX, $actualY-$dasSUMY"
		#read -n1 kbd

		echo ">nothing" > tempfastas/Y_${counterY}.fasta	
		(tail -c +"$actualY" "$seqNameY.fix" | head -c "$ratioY") >> tempfastas/Y_${counterY}.fasta
	
		echo ">nothing" > tempfastas/Y_${counterY}_rev.fasta	
		(tail -c +"$fakeY" "$seqNameY.fix" | head -c "$ratioY") >> tempfastas/Y_${counterY}_rev.fasta

		counterY=`expr $counterY + 1`
		
	fi
	
	
	
	# switcher
	if [[ turn -eq 0 ]]; then
		turn=1
	else
		turn=0
		
		counterXprev=`expr $counterX - 1`
		counterYprev=`expr $counterY - 1`
		
		# run gecko on this comparison



		#echo "($BINDIR/gecko tempfastas/X_${counterXprev}.fasta tempfastas/Y_${counterYprev}.fasta temp.frags $LEN $SIM $WL)"
		($BINDIR/gecko tempfastas/X_${counterXprev}.fasta tempfastas/Y_${counterYprev}.fasta temp.frags $LEN $SIM $WL f) | grep "Frags in"
		#echo "$BINDIR/filterFrags temp.frags $LEN $SIM > X_${counterXprev}-Y_${counterYprev}.csv"
		($BINDIR/filterFrags temp.frags $LEN $SIM > X_${counterXprev}-Y_${counterYprev}.csv)
		#cat X_${counterXprev}-Y_${counterYprev}.csv
		#read -n1 kbd
		($BINDIR/gecko tempfastas/X_${counterXprev}.fasta tempfastas/Y_${counterYprev}_rev.fasta temp.frags $LEN $SIM $WL r)  | grep "Frags in"
		($BINDIR/filterFrags temp.frags $LEN $SIM > X_${counterXprev}-Y_${counterYprev}_rev.csv)

		#(tail -n +18 X_${counterXprev}-Y_${counterYprev}.csv | awk -F "," -v OFS=',' -v a="$actualX" -v b="$actualY" '{print $1,$2+a,$3+b,$4+a,$5+b,$6,$7,$8,$9,$10,$11,$12,$13,$14}') >> all-results/master.csv


		(tail -n +18 X_${counterXprev}-Y_${counterYprev}.csv | awk -F "," -v OFS=',' -v a="$actualX" -v b="$actualY" '{if($6 == "f") print $1,$2+a,$3+b,$4+a,$5+b,$6,$7,$8,$9,$10,$11,$12,$13,$14; else print $1,$2+a,$5+b,$4+a,$3+b,$6,$7,$8,$9,$10,$11,$12,$13,$14 ;}') >> all-results/master.csv
		(tail -n +18 X_${counterXprev}-Y_${counterYprev}_rev.csv | awk -F "," -v OFS=',' -v a="$actualX" -v b="$fakeY" -v ratio="$ratioY" '{if($6 == "f") print $1,$2+a,(ratio-$3)+b,$4+a,(ratio-$5)+b,$6,$7,$8,$9,$10,$11,$12,$13,$14; else print $1,$2+a,(ratio-$3)+b,$4+a,(ratio-$5)+b,$6,$7,$8,$9,$10,$11,$12,$13,$14 ;}') >> all-results/master.csv

		


		rm -rf temp.frags X_${counterXprev}-Y_${counterYprev}.csv X_${counterXprev}-Y_${counterYprev}_rev.csv



		
	fi
	
done
rm $seqNameX.fix $seqNameY.fix


