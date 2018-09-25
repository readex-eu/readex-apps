find . -type f -iname "*.csv" -print0 | while IFS= read -r -d $'\0' line 
do 
echo "$line"; 
echo 'Energy----------------------'
sed -n '/# Blade summary/{n;p}' $line | sed -n 's/Energy consumption \[J\] (Samples),\(.*\)/\1/p'
#echo 'Runtime----------------------'
#sed -n '/# Job info - hdeem/{N;N;N;p}' $line | sed -n 's/Runtime of function \[s\],\(.*\)/\1/p'
done

find . -type f -iname "*.csv" -print0 | while IFS= read -r -d $'\0' line 
do 
echo "$line"; 
#echo 'Energy----------------------'
#sed -n '/# Blade summary/{n;p}' $line | sed -n 's/Energy consumption \[J\] (Samples),\(.*\)/\1/p'
echo 'Runtime----------------------'
sed -n '/# Job info - hdeem/{N;N;N;p}' $line | sed -n 's/Runtime of function \[s\],\(.*\)/\1/p'
done

