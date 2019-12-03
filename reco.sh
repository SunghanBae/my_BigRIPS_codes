#!/bin/bash
ncpu=4

myfunc(){
(
cd $BigDIR
date +"%d, %T"

nam=run${1}.root
if [ -e $nam ]; then
		echo "/reco ${1} ${2}"
		./reco ${1} ${2} >> ./log/reco_${1}_tof_${2}
		echo $nam 'PID reconstruction done'
else echo "$nam does not exist"
fi
)
}

export -f myfunc

#for ((i=1185;i<=1263;i+=1))
#for i in {1185..1201} {1252..1263}
for i in {1248..1250}
do
#for 1252~1263:      445.18
#for 119x, n-rich:   445.00
if [ ${i} -lt 1248 -a ${i} -ge 1184 ]; then
		toff=445.00
		else if [ ${i} -ge 1248 -a ${i} -lt 1264 ]; then
		toff=445.18
		else if [ ${i} -ge 1181 -a ${i} -lt 1184 ]; then
		toff=447.5283
		else echo "The run number is weird!!" 
			break
fi
fi
fi

list=${list}"${i} ${toff}\n"
#echo ${i}" "${toff}

done

echo -e $list | xargs -P$ncpu -I{} bash -c "myfunc {}"
