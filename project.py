# -*- coding: utf-8 -*-
"""
Created on Sun May  3 23:47:35 2015

@author: Asus
"""
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
import re
import mysql.connector

def main():
    file = open("blast.txt", "w").close()
    file = open("blast.txt", "a")
    bestand = "seqschamppro.txt"    #tab delimited text file van excel bestand 100 sequenties groep 1
    b = open(bestand, "r")
    a = ""
    count = 1
    seqs = []
    for regel in b:
        a+=regel
    i = a.split("@HWI")  
    i.remove("")
    
    d1, d2 = dikt(i)

          
    for d1key, d1value in d1.items():
    	results = []
	    bool = False
        if count < 102:
            blastnr =  ("\n\n\n------------blast number: "+str(count)+"-----------------------------------------------\n\n\n")
            idheader =  ("HEAD:      "+d1key)
            seqheader = ("SEQUENCE:  "+d1value+"\n\n")
            print (blastnr)
            print (idheader)
            print (seqheader)			
	        file.write("\n"+str(blastnr))
	        file.write("\n"+str(idheader))
	        file.write("\n"+str(seqheader))
            source_seq(bool, d1key, d1value)
            result = blast(d1value, file, d1key)
            results.append(result)
            boolean(result, d1key)
            if result == False:
                noresult =  ("no results found for blast number "+str(count))
        		print (noresult)
        		count += 1
            else:
                endresults =  ("\n\n\n------------end of results blast number "+str(count)+"-------------------------------------")
		        print (endresults)
		        file.write("\n"+str(endresults))
                count += 1
		    bool = True
	
        else:
            print ("done")

def boolean(b, head):
     if b == False: 
        a = 0
     elif b == True: 
        a = 1
     conn = mysql.connector.connect(host="ithurtswhenip.nl",user = "richard", password = "richard", db = "blast", port = 3307)
     cursor = conn.cursor()
     query = ("UPDATE `blast`.`SOURCE_SEQ` SET `match-found-bool`='"+str(a)+"' WHERE `seq-id`='"+str(head)+"';")   
     print (query)
     cursor.execute(query)
     conn.commit()
     cursor.close()
     conn.close()     

def source_seq(bool, head, seq): 
     conn = mysql.connector.connect(host="ithurtswhenip.nl",user = "richard", password = "richard", db = "blast", port = 3307)
     cursor = conn.cursor()
     query = ("INSERT INTO `blast`.`SOURCE_SEQ`(`seq-id`, `match-found-bool`, `nucleotide-sequence`) VALUES ('"+str(head)+"', "+str(bool)+",'"+str(seq)+"');")
     cursor.execute(query)   
     conn.commit()
     print (query)
     cursor.close()
     conn.close()

def blast(seq, file, head):
    bool = True	
    results_handle = NCBIWWW.qblast("blastx", "nr", seq, hitlist_size=10, expect=1)
    blast_results = results_handle

    blast_records = NCBIXML.parse(results_handle)

    E_VALUE_THRESH = 1

    c = 1
    for blast_record in blast_records:
        desc = blast_record.descriptions
        if desc == []:
    	    bool = False
            else: 
            numhits = len(desc)
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                if hsp.expect < E_VALUE_THRESH:
                    alignmentnr =  ("-----Alignment-"+str(c)+"----")
                    hithead =  alignment.title
                    length = alignment.length
                    evalue =  hsp.expect
                    gaps =  hsp.gaps
                    query =  hsp.query
                    sub = hsp.sbjct
                    match = hsp.match
                    score = hsp.score
                    identity = hsp.identities
    
                    querylen = len(query)                 
                    qcov = float(identity) / float(querylen) * float(100)
		            if numhits == None:
		                numhits ="NULL"
                    print ("hithead:   ", hithead)
                    print ("length:    ", length)
                    print ("e-value:   ", evalue)
                    print ("gaps:      ", gaps)
                    print (query)
                    print (match)
                    print (sub)
                    print ("identity:  ", identity)
                    print ("score:     ", score)
		            print ("qcov:      ", qcov)
		            print (len(query))
                    # r =  ("sequence:  ", sequence,"\nlength:    ", length,"\ne-value:   ", evalue,"\ngaps:      ", gaps, "\n",query, "\n",match, "\n",sub, "\nidentity:  ", identity, "\nscore:     ", score)
                    # print (r)
                    file.write("\n"+str(alignmentnr))
                    file.write("\n\n HITHEAD:  \n"+str(hithead))
    	            file.write("\n LENGTH:         "+str(length))
    	            file.write("\n E-VALUE:        "+str(evalue))
    	            file.write("\n GAPS:           "+str(gaps))
    	            file.write("\n "+str(query))
    	            file.write("\n "+str(match))
    	            file.write("\n "+str(sub))
    	            file.write("\n IDENTITY%:      "+str(identity))
    	            file.write("\n SCORE:          "+str(score))
    	            file.write("\n NUMBER OF HITS: "+str(numhits)+"\n")
#		        print ("\n***********************************************\n")
#		        print ("head\n"+str(head))
#		        print ("\n***********************************************\n")
#		        print ("numhits\n"+str(numhits))
#		        print ("\n***********************************************\n")
#		        print ("sequence\n"+str(hithead))
#		        print ("\n***********************************************\n")
#		        print ("score\n"+str(score))
#		        print ("\n***********************************************\n")
#		        print ("qcov\n"+str(qcov)+" start: "+str(querys)+" end: "+str(querye)+" len: "+str(len(query)))
#		        print ("\n***********************************************\n")
#		        print ("evalue\n"+str(evalue))
#		        print ("\n***********************************************\n")
#		        print ("identity\n"+str(identity))
# 		        print ("\n***********************************************\n")
	            if c == 1:
		           # print ("INSERT INTO `blast`.`BLAST_RESULT`(`amount-hits`, `result1-accession-code`, `result1-alignment-score`, `result1-query-coverage-percentage`, `result1-E-value`, `result1-ident-percentage`, `SOURCE_SEQ_seq-id`) VALUES ("+str(numhits)+", '"+str(hithead)+"', "+str(score)+", "+str(qcov)+", "+str(evalue)+", "+str(identity)+",'"+str(head)+"');")
	                conn = mysql.connector.connect(host="ithurtswhenip.nl",user = "richard", password = "richard", db = "blast", port = 3307)
		            cursor = conn.cursor()
		            query1 = ("INSERT INTO `blast`.`BLAST_RESULT`(`amount-hits`, `result1-accession-code`, `result1-alignment-score`, `result1-query-coverage-percentage`, `result1-E-value`, `result1-ident-percentage`, `SOURCE_SEQ_seq-id`) VALUES ("+str(numhits)+", '"+str(hithead)+"', "+str(score)+", "+str(qcov)+", "+str(evalue)+", "+str(identity)+",'"+str(head)+"');")
		           # rows = str(cursor.fetchall())	
	               #  print (rows)
		            cursor.execute(query1)
		            print query1
                    conn.commit()
               	    cursor.close()
		            conn.close()
 		            c += 1
		        else: 
                    conn = mysql.connector.connect(host="ithurtswhenip.nl",user = "richard", password = "richard", db = "blast", port = 3307)
                    cursor = conn.cursor()
                    query2 = ("UPDATE `blast`.`BLAST_RESULT` SET `result"+str(c)+"-accession-code`='"+str(hithead)+"', `result"+str(c)+"-alignment-score`='"+str(score)+"', `result"+str(c)+"-query-coverage-percentage`='"+str(qcov)+"', `result"+str(c)+"-E-value`='"+str(evalue)+"', `result"+str(c)+"-ident-percentage`='"+str(identity)+"' WHERE `SOURCE_SEQ_seq-id`='"+str(head)+"';")
                   # query2 = ("INSERT INTO `blast`.`BLAST_RESULT`(`result"+str(c)+"-accession-code`, `result"+str(c)+"-alignment-score`, `result"+str(c)+"-query-coverage-percentage`, `result"+str(c)+"-E-value`, `result"+str(c)+"-ident-percentage`) VALUES ('"+str(hithead)+"', "+str(score)+", "+str(qcov)+", "+str(evalue)+", "+str(identity)+");")               
		            cursor.execute(query2)
                    print (query2)
                    conn.commit()
                    cursor.close()
                    conn.close()
		            c += 1
    return (bool)
def dikt(i):
    dikt1 = {}
    dikt2 = {}
    for item in i:     
        lijst = item.split("\t")
        head = lijst[0]
        
        p1 = re.compile("_1")
        m1 = p1.search(head)
        p2 = re.compile("_2")
        m2 = p2.search(head)
        if m1:
            head1 = ("@HWI"+lijst[0])
            seq1 = lijst[1]
            score1 = lijst[2]            
            dikt1.update({head1:seq1})
        if m2: 
            head2 = ("@HWI"+lijst[0])
            seq2 = lijst[1]
            score2= lijst[2]            
            dikt2.update({head2:seq2})
        
    return(dikt1, dikt2)
	
# def dikt(i):
    # dikt = {}
   
    # for item in i:     
        # lijst = item.split("\t")
                
        # head = ("@HWI"+lijst[0])
        # seq = lijst[1]
        # score = lijst[2]            
        # dikt.update({head:seq})
      
        
    # return(dikt)

# def newfile(sequence, length, evalue, gaps, query, match, sub, identity, score):
    # count = 0
    # file = open("blast.txt", "a")
    # if count = 0: 
		# file.write("RESULTS\n")
    # else:    
		# file.write("\n SEQUENCE \n "+str(sequence))
		# file.write("\n LENGTH:     "+str(length))
		# file.write("\n E-VALUE:    "+str(evalue))
		# file.write("\n GAPS:       "+str(gaps))
		# file.write("\n"+str(query))
		# file.write("\n"+str(match))
		# file.write("\n"+str(sub))
		# file.write("\n IDENTITY%:  "+str(identity))
		# file.write("\n SCORE:      "+str(score)+"\n")

# def update_progress(progress):
#	print ('\r[{0}] {1}%'.format('#'*(progress/10), progress)

main()
