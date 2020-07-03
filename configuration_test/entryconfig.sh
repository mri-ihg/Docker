# Prepares the XML config files

CONFIG=fields;
CONFIGPATH="."

# Create target files
TARGETFILES="";
for TEMPLATE in $CONFIGPATH/*template*.xml
do
	TARGETFILE=$( echo ${TEMPLATE} | sed s/\.template//g);
	TARGETFILES=${TARGETFILE}" "${TARGETFILES}
	cp ${TEMPLATE} ${TARGETFILE};
done


# Write KEY/VAL on config files
while read ENTRY
do
	KEY=$( echo $ENTRY | awk -F"=" '{print $1}' | awk '{print $1}'  | tr -d " " );
	VAL=$( echo $ENTRY | awk -F"=" '{print $2}' | awk '{print $1}' | awk -F";" '{print $1}' | tr -d " " );

	echo "$KEY : $VAL;";

	for TARGETFILE in ${TARGETFILES}
	do
		sed -i s/${KEY}/${VAL}/g ${TARGETFILE}
	done

done < <(cat $CONFIG | grep "="  |  awk -F"#" '{print $1}' | grep -ve "^$"; );
