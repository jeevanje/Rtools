Rfile=$1
cp $Rfile ${Rfile}_ncdf
sed -e 's|get.var.ncdf|ncvar_get|'\
    < ${Rfile}_ncdf > ${Rfile}
rm ${Rfile}_ncdf
