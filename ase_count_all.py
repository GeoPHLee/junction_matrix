#A3SS
#chr18:13114129:13114250:+@chr18:13116373|13116376:13116502:+
#chr3:7620109:7621044:+@chr3:7721712|7721736:7721982:+
#chr3:184019325:184019446:+@chr3:184019559|184019646:184019859:+
#chr18:693882:693908:-@chr18:690631|690745:690549:-
#SE
#chr18:2544652:2544758:-@chr18:2544194:2544285:-@chr18:2537524:2539144:-


#read_list
read_list1=[(18,13116373,13116502),(18,13116373,13116502),(3,7721712,7721736),(3,7721712,7721736),(18,13116373,13116502),(18,13116373,13116502),(18,13116373,13116502)]
read_list2=[(18,13116373,13116502),(3,184019559,184019859),(3,184019559,184019646),(3,184019559,184019646),(3,184019559,184019859),(18,13116373,13116502),(18,13116373,13116502),(18,13116373,13116502),(3,184019559,184019859),(18,13116373,13116502),(3,7721712,7721982),(3,7721712,7721982),(3,7721712,7721982)]
read_list3=[(18,13116373,13116502),(3,184019559,184019646),(3,184019559,184019646),(3,184019559,184019859),(3,184019559,184019646),(3,184019559,184019647),(3,184019559,184019646),(3,184019559,184019859),(18,13116373,13116502),(18,13116373,13116502),(3,7721712,7721736),(3,7721712,7721536),(18,13116373,13116502),(3,184019559,184019859),(18,13116373,13116502),(3,7721712,7721982),(3,7721712,7721982),(3,7721712,7721982)]

ase_list={'ase1-in':(18,13116373,13116502),
'ase1-out':(18,13116373,13116376),
'ase2-in':(3,7721712,7721736),
'ase2-out':(3,7721712,7721982),
'ase3-in':(3,184019559,184019646),
'ase3-out':(3,184019559,184019859),
}


ase_list={
(18,13116373,13116502):("ase1","in"),
(18,13116373,13116376):("ase1","out"),
(3,7721712,7721736):("ase2","in"),
(3,7721712,7721982):("ase2","out"),
(3,184019559,184019646):("ase3","in"),
(3,184019559,184019859):("ase3","out"),
(18,690631,690549):'ase4-in',
(18,690745,690549):'ase4-out',
}

ase_cout={
"ase1":{"in":count,"out":count},
"ase2":{"in":count,"out":count},
"ase3":{"in":count,"out":count},
"ase4":{"in":count,"out":count}
}

#ase_list
#(chrosome,start,end,strand):("ase_name",in_or_out)
#A3SS/A5SS in=anchor1 out=anchor2
#SE in=anchor1 out={anchor}


junction=(chrosome,start,end,strand)
junction:("ase_name",in_or_out)

#ase_cout
#"ase_name":{"in":count,"out":count}
in_out_counter={"in":count,"out":count}
ase_name:in_out_counter

#read_count
#ase_cout[ase_list[parsed_read][0]][ase_list[parsed_read][1]]=ase_cout[ase_list[parsed_read][0]][ase_list[parsed_read][1]]+1

for parsed_read in read_list3:
	#query_result=("ase1","in"or"out) if the ase in the table,else raise key error
	try:
		query_result=ase_list[parsed_read]
		ase_name=query_result[0]
		junction_anchor=query_result[1]
		ase_cout[ase_name][junction_anchor]=ase_cout[ase_name][junction_anchor]+1
	except KeyError as e:
		print("no this junction")
#		print(parsed_read)


#compute psi
[x["in"]/(x["out"]+x["in"]) for psi_counter in ase_cout.values() if x["out"]+x["in"]!=0 ]

#Lists_of_ASE
ase_list={
	
}

junctions_list=[]
for read in junction_list:
           print(read.cigar,read.pos)
           read_start=read.pos
           junction_flag=0
           #read.cigar[0][1],read.pos+read.cigar[0][1]+read.cigar[1][1]
           for part in read.cigar:
                part_status=part[0]
                part_length=part[1]
                if part_status==0:
                                junction_start=read_start+part_length
                                junction_flag=1
                                #print(junction_start)
                if part_status==3:
                                junction_end=junction_start+part_length
                                junction=(junction_start,junction_end)
                                junctions_list.append((junction_start,junction_end))
                                junction_flag=0

import gffutils ,Bio
db = gffutils.create_db("/data-2/reference/hg19/annotation/gencode.v28lift37.annotation.gtf", "hg19")
print("OK")


chr7:73604152:73604248:+@chr7:73604577:73604636:+@chr7:73609071:73609208
samfile = pysam.AlignmentFile("CHG026340-WYB-1-S1-ATCACG_L002/alignments.sort.bam", "rb")

 read_get=samfile.fetch('chr7', 73604002, 73609358,"+")

 junction_list=[read for read in read_get if read.cigar!=[(0,100)]]

chr7:73604152:73604248:+@chr7:73604577:73604636:+@chr7:73609071:73609208:+ EIF4H 0 0 0.398956 nan