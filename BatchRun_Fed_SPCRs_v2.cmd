cd Dropbox (Gladstone)/eVOLVER_and_Retro_Record/MiSeq_Data/
for /L %%i in (1,1,9) do python Fed_SPCRs_v2.py msSBK_2 msSBK_2_0%%i_trimmed local
for /L %%i in (10,1,42) do python Fed_SPCRs_v2.py msSBK_2 msSBK_2_%%i_trimmed local

cmd /k

REF python Fed_SPCRs_v2.py msSBK_1 msSBK_1_56_trimmed local
