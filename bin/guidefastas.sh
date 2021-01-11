#!/bin/bash 

BINDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

if [ $# -ne 9 ]; then
   echo " ==== ERROR ... you called this script inappropriately."
   echo ""
   echo "   usage:  $0 seqXName.fasta seqYName.fasta guideFile dimension LEN SIM WL print[0/1] ids/names[0/1]"
   echo ""
   exit -1
fi

binary_search(){
  TARGET=$1
	name=$2[@]
	TO_SEARCH=("${!name}")
	LENGTH=${#TO_SEARCH[@]}

	START=0
	END=$((LENGTH - 1))
	while [[ $START -le $END ]]; do
		MIDDLE=$((START + ((END - START)/2)))
		ITEM_AT_MIDDLE=${TO_SEARCH[MIDDLE]}
		#echo "comparing $TARGET with $ITEM_AT_MIDDLE"
		if [[  $ITEM_AT_MIDDLE -gt $TARGET ]]; then
			END=$((END-MIDDLE-1))
		elif [[ $ITEM_AT_MIDDLE -lt $TARGET ]]; then
			START=$((MIDDLE+1))
		else
			echo $MIDDLE
			return 0
		fi
	done
	if [[ $TARGET -gt ITEM_AT_MIDDLE ]]; then
		echo $MIDDLE
	else
		echo $((MIDDLE-1))
	fi
	return 0
}




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

### Build index to account on which sequence we are working on ################

#mapfile -t indexX < <(grep ">" $1 -b | awk -F ":" '{print $1}')
#mapfile -t indexY < <(grep ">" $2 -b | awk -F ":" '{print $1}')
grep ">" $1 -b | awk -F ":" '{first=$1; $1=""; print first $0;}' | sed 's/ /_/g' | sed 's/>//g' > $seqNameX.indices
grep ">" $2 -b | awk -F ":" '{first=$1; $1=""; print first $0;}' | sed 's/ /_/g' | sed 's/>//g' > $seqNameY.indices
#grep ">" $1 -b | awk -F ":" '{print $1 $2}' > $seqNameX.indices
#grep ">" $2 -b | awk -F ":" '{print $1 $2}' > $seqNameY.indices

#closest=$(binary_search 649261907 indexX)

##################

guided=$3
dimension=$4

LEN=$5
SIM=$6
WL=$7
print=$8
nameIDs=$9

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
	
		actualX=`expr $i - 0` 
		actualX=$(($actualX * $ratioX))

	
		echo ">nothing" > tempfastas/X_${counterX}.fasta
		(tail -c +"$actualX" "$seqNameX.fix" | head -c "$ratioX") >> tempfastas/X_${counterX}.fasta
		
	
		counterX=`expr $counterX + 1`
	fi
	
	if [[ turn -eq 1 ]]; then
		
		actualY=`expr $i - 0` # used to be -2
		actualY=$(($actualY * $ratioY))
		fakeActualY=`expr $i + 0`
		fakeY=$(($fakeActualY * $ratioY))
		
		dasSUMX=`expr $actualX + $ratioX`
		dasSUMY=`expr $actualY + $ratioY`
		echo "(i: $i) coord X,Y: $actualX-$dasSUMX, $actualY-$dasSUMY"
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



		($BINDIR/gecko tempfastas/X_${counterXprev}.fasta tempfastas/Y_${counterYprev}.fasta temp.frags $LEN $SIM $WL f) | grep "Frags in"
		($BINDIR/filterFrags temp.frags $LEN $SIM > X_${counterXprev}-Y_${counterYprev}.csv)
		
		($BINDIR/gecko tempfastas/X_${counterXprev}.fasta tempfastas/Y_${counterYprev}_rev.fasta temp.frags $LEN $SIM $WL r)  | grep "Frags in"
		($BINDIR/filterFrags temp.frags $LEN $SIM > X_${counterXprev}-Y_${counterYprev}_rev.csv)


		#echo "$BINDIR/gecko tempfastas/X_${counterXprev}.fasta tempfastas/Y_${counterYprev}_rev.fasta temp.frags $LEN $SIM $WL r"
		#echo "$BINDIR/filterFrags temp.frags $LEN $SIM > X_${counterXprev}-Y_${counterYprev}_rev.csv"



		(tail -n +18 X_${counterXprev}-Y_${counterYprev}.csv | awk -F "," -v OFS=',' -v a="$actualX" -v b="$actualY" '{if($6 == "f") print $1,$2+a,$3+b,$4+a,$5+b,$6,$7,$8,$9,$10,$11,$12,$13,$14; else print $1,$2+a,$5+b,$4+a,$3+b,$6,$7,$8,$9,$10,$11,$12,$13,$14 ;}') >> all-results/master.csv

		# Extract forward frags

		for thisfrag in $(tail -n +18 X_${counterXprev}-Y_${counterYprev}.csv | awk -F "," -v OFS=',' -v a="$actualX" -v b="$actualY" '{if($6 == "f") print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14; else print $1,$2,$5,$4,$3,$6,$7,$8,$9,$10,$11,$12,$13,$14 ;}')
		do
			#echo "this line $thisfrag"
			arrIN=(${thisfrag//,/ })

			xStart=${arrIN[1]}
			yStart=${arrIN[2]}
			xEnd=${arrIN[3]}
			yEnd=${arrIN[4]}

			# Locate in which sequence they are
			#seqXid=$(binary_search $xStart indexX)
			#seqYid=$(binary_search $yStart indexY)
			#echo "Located $xStart at $seqXid with value ${indexX[$seqXid]}"

			#echo $thisfrag
			# Write the fragments
			#(echo $thisfrag | awk -F "," -v OFS=',' -v a="$actualX" -v b="$actualY" -v xname="$seqXid" -v yname="$seqYid" '{if($6 == "f") print $1,$2+a,$3+b,$4+a,$5+b,$6,$7,$8,$9,$10,$11,$12,xname,yname; else print $1,$2+a,$5+b,$4+a,$3+b,$6,$7,$8,$9,$10,$11,$12,xname,yname ;}') >> all-results/master.csv

			if [[ $print -eq 1 ]]; then

				awk -v x1="$xStart" -v x2="$xEnd" -v y1="$yStart" -v y2="$yEnd" 'NR==FNR && FNR==2{l1=$0} NR!=FNR && FNR==2{l2=$0} 
				BEGIN{ x1++; y1++; x2++; y2++; }
	        	       	END{good=0; split(l1, s1, ""); split(l2, s2, "");
				for(i=x1; i<=x2; i++){ printf("%c", s1[i]); } printf("\n");
				j = y1;
		                for(i=x1; i<=x2; i++){ if(s1[i]==s2[j]) {good++; printf("|");} else printf(" "); j++;  }; printf("\n");
				for(i=y1; i<=y2; i++){ printf("%c", s2[i]); } printf("\n");
       	        		print "@ FORWARD STRAND Identity:",good "/" x2-x1+1, "("100*good/(x2-x1+1)"%)";  }' tempfastas/X_${counterXprev}.fasta tempfastas/Y_${counterYprev}.fasta
			fi

		done

		
		(tail -n +18 X_${counterXprev}-Y_${counterYprev}_rev.csv | awk -F "," -v OFS=',' -v a="$actualX" -v b="$fakeY" -v ratio="$ratioY" '{if($6 == "f") print $1,$2+a,(ratio-$3)+b,$4+a,(ratio-$5)+b,$6,$7,$8,$9,$10,$11,$12,$13,$14; else print $1,$2+a,(ratio-$3)+b,$4+a,(ratio-$5)+b,$6,$7,$8,$9,$10,$11,$12,$13,$14 ;}') >> all-results/master.csv

		# Extract reverse frags

		for thisfrag in $(tail -n +18 X_${counterXprev}-Y_${counterYprev}_rev.csv | awk -F "," -v OFS=',' -v a="$actualX" -v b="$fakeY" -v ratio="$ratioY" '{if($6 == "f") print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14; else print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14 ;}')
		do

			#echo "this line $thisfrag"
			arrIN=(${thisfrag//,/ })
			xStart=${arrIN[1]}
			yStart=${arrIN[2]}
			xEnd=${arrIN[3]}
			yEnd=${arrIN[4]}

			# Locate in which sequence they are
            #seqXid=$(binary_search $xStart indexX)
            #seqYid=$(binary_search $yStart indexY)
            #echo "Located $xStart at $seqXid with value ${indexX[$seqXid]}"

            #echo $thisfrag
            # Write the fragments
            #(echo $thisfrag | awk -F "," -v OFS=',' -v a="$actualX" -v b="$fakeY" -v ratio="$ratioY" -v xname="$seqXid" -v yname="$seqYid" '{if($6 == "f") print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,xname,yname; else print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,xname,yname ;}') >> all-results/master.csv


			if [[ $print -eq 1 ]]; then

				awk -v x1="$xStart" -v x2="$xEnd" -v y1="$yStart" -v y2="$yEnd" -v ytotal="$ratioY" 'NR==FNR && FNR==2{l1=$0} NR!=FNR && FNR==2{l2=$0}
				BEGIN{ x1++; x2++; y1 = ytotal - y1; y2 = ytotal - y2; y1--; y2--; print x1,x2,y1,y2 }
				END{good=0; split(l1, s1, ""); split(l2, s2, "");
				for(i=x1; i<=x2; i++){ printf("%c", s1[i]); } printf("\n");

				j = y1;
				for(i=x1; i<=x2; i++){ 
					charS2 = "N";
					if(s2[j] == "A") charS2 = "T";
					if(s2[j] == "C") charS2 = "G";
					if(s2[j] == "G") charS2 = "C";
					if(s2[j] == "T") charS2 = "A";

					if(s1[i] == charS2) {good++; printf("|");} else printf(" "); 
					j--;  
				}; printf("\n");

				for(i=y1; i>=y2; i--){
					if(s2[i] == "A") printf("%c", "T"); 
					if(s2[i] == "C") printf("%c", "G"); 
					if(s2[i] == "G") printf("%c", "C"); 
					if(s2[i] == "T") printf("%c", "A"); 
					if(s2[i] == "N") printf("%c", "N"); 
				} printf("\n");
				print "@ REVERSE STRAND Identity:",good "/" x2-x1+1, "("100*good/(x2-x1+1)"%)";  }' tempfastas/X_${counterXprev}.fasta tempfastas/Y_${counterYprev}_rev.fasta
			fi

		done

		


		rm -rf temp.frags X_${counterXprev}-Y_${counterYprev}.csv X_${counterXprev}-Y_${counterYprev}_rev.csv



		
	fi
	
done
rm $seqNameX.fix $seqNameY.fix

lx=$(wc -l ${seqNameX}.indices)
ly=$(wc -l ${seqNameY}.indices)

$BINDIR/masterToNames all-results/master.csv $seqNameX.indices $seqNameY.indices $lx $ly $nameIDs

#rm $seqNameX.indices $seqNameY.indices


