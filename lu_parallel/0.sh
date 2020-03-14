path=$1
files=$(ls $path)

for filename in $files
do
	./demos $path/$filename 
done
