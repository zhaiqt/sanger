/nhome/zhaiqi1/NGS/IgBlast/1.4.0/ncbi-igblast-1.4.0/bin/igblastp -germline_db_V /home/zhaiqi1/NGS/IgBlast/database/human_gl_V  -domain_system imgt -num_alignments_V 1 -query assembled_protein  -organism human > IgBlastp.txt 

../../IgBlast/1.4.0/ncbi-igblast-1.4.0/bin/igblastp -germline_db_V /NGS/IgBlast/database/human_gl_V  -domain_system imgt -num_alignments_V 1 -query assembled_protein  -organism human > IgBlastp.txt



/home/zhaiqi1/NGS/IgBlast/igblastp -germline_db_V /NGS/IgBlast/database/human_gl_V  -domain_system imgt -num_alignments_V 1 -query assembled_protein  -organism human > IgBlastp.txt


python /home/zhaiqi1/NGS/IgBlastp_wrapper.py -s human -i /home/zhaiqi1/NGS/RESULTS/protein.txt


python /home/zhaiqi1/NGS/IgBlastp_wrapper.py -s human -i /home/zhaiqi1/NGS/RESULTS/assembled_protein  
python WrapIgBlastn.py  -s mouse -i test/tmp


qsub -b y python /home/zhaiqi1/NGS/mycode/Ab_NGS_4/PostAnalysis.py  -d /dlab/NGS/usem-seqanalysis/160808_zhaiqi1_miseq_CDH17/H-easy/RESULTS -c H  -s mouse




qsub -b y python PostAnalysis.py  -d /dlab/NGS/usem-seqanalysis/160314_zhaiqi1_miseq_HBx52-60DNA.20160214_AN2N4/LC/RESULTS -c L  -s mouse

