for salt in 0.*
do
	cd $salt
	pwd
	python ../density.py cores.dump &
	cd -
done
wait
