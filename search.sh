# inputs:
# 	genome to search for frameshifts
# 	db file
# 	Anti SD sequence
# 	number of upstream bases to consider for stemloop/ pseudoknot

#declare default values
upstream=14 #number of upstream bases to evaluate for arrest sequence
downstream=81 #number of downstream bases to evaluate for conserved structure
sd_frame=15 #number of upstream bases to consider for SD-like sequence
structure_start=11 #number of bases after beginning of slippery sequence to start looking for structure

while getopts o:q:d:l:u:s:r:p:i: flag
		do
		        case "${flag}" in
		                o) out=${OPTARG};;
		                q) fa=${OPTARG};;
		                d) db=${OPTARG};;
						l) downstream=${OPTARG};;
						u) upstream=${OPTARG};;
						s) sd_frame=${OPTARG};;
						r) rscape_path=${OPTARG};;
						p) path_to_prf_search=${OPTARG};;
						i) structure_start=${OPTARG};;
		                *) echo "Invalid option: -$flag" ;;
		        esac
		done

 
# up=$(($upstream % 3))
# if [[ $up -ne 0 ]] ; then 
# 	echo "Error: Upstream bases (-u) must be divisible by 3."
# 	exit 1
# fi

time_stamp=$(date "+%Y_%m_%d-%H_%M_%S")
root=$out$time_stamp

python3 split_fa.py -F ${fa} -O ${root}
wait

python3 calculate_model_params.py -F /Users/mollysacks/thesis/prf-search/All-genes-of-E.-coli-K-12-substr.-MG1655.fasta
wait

#count number of directories created
total=`find ${root} -maxdepth 1 -type d | wc -l | xargs`

python3 slippery_sequence_search.py -D ${root} -L ${downstream} -U ${upstream} -B ${structure_start}
wait

find_features () {
	for loc in ${cDNA}/*
	# for each slippery sequence
	do
		if [[ -d $loc ]]
        then
        	# find fasta file
        	for seq in $loc/*.fa
        	do
        		# build multiple alignment
        		nhmmer -A ${seq}.sto -o ${seq}.txt ${seq} ${db}
        		
        		# make sure there were hits before proceeding
        		if [ -s ${seq}.sto ]
        		then
        			# continue building alignment
        			hmmbuild --fast --symfrac 0.0 -o ${seq}.iter1.hmm.txt ${seq}.iter1.hmm ${seq}.sto
					nhmmer -A ${seq}.iter1.sto -o ${seq}.iter1.txt ${seq}.iter1.hmm ${db}
					hmmbuild --fast --symfrac 0.0 -o ${seq}.iter2.hmm.txt ${seq}.iter2.hmm ${seq}.iter1.sto
					nhmmer -A ${seq}.iter2.sto -o ${seq}.iter2.txt ${seq}.iter2.hmm ${db}
					hmmbuild --fast --symfrac 0.0 -o ${seq}.iter3.hmm.txt ${seq}.iter3.hmm ${seq}.iter2.sto
					nhmmer -A ${seq}.iter3.sto -o ${seq}.iter3.txt ${seq}.iter3.hmm ${db}
					
					t=`python3 mask_params.py -F ${seq} -S ${seq}.iter3.sto`
					ti=$(($t +26 ))

					esl-alimask -t -o ${seq}.final.sto ${seq}.iter3.sto ${t}: > ${seq}.final.txt
					esl-alimask -t -o ${seq}.structure.sto ${seq}.iter3.sto ${ti}: > ${seq}.structure.txt
					esl-alistat --cinfo ${seq}.conservation ${seq}.final.sto > ${seq}.conservation.txt 
					
					# fold
					mkdir ${loc}/CaCoFold
					cd ${loc}/CaCoFold
					${rscape_path} --nofigures --fold ${path_to_prf_search}/${seq}.structure.sto &> stdout
					cp *.fold.sto ${path_to_prf_search}/${loc}
					cd ..

					# check for structure and look for other (optional) features
					python3 ${path_to_prf_search}/generate_loc_report.py \
						-F *.fa \
						-S *.fold.sto \
						-B 5 \
						-R ${sd_frame} \
						-P 5 \
						-O ${path_to_prf_search}/${root}.report \
						-L ${loc} \
						-M ${structure_start} \
						-Q ${path_to_prf_search}/frequencies \
						-C *.conservation

					# Get conservation info

					cd ${path_to_prf_search}
				else
        			#echo "No hits satisfy inclusion thresholds; no alignment saved for ${loc}"
        			#echo "Removing ${loc} from candidate PRF sites"
        			rm -R $loc
        			continue
     			fi
				done
        fi
	done
}

i="1"
sss="0"

for cDNA in ${root}/*
do 
	j=$(($i % 5))
	if [[ $j -ne 0 ]] ; 
	then 
		find_features ${cDNA} &
	else
		find_features ${cDNA}
		wait
		python3 reformat_output.py -O ${root}.report
	fi
	i=$(( i + 1 ))
done
wait

python3 random_probability_eval.py -I ${fa}

python3 create_plots.py -O ${root}.report.tsv -N null.tsv -C ${fa}.results.tsv




