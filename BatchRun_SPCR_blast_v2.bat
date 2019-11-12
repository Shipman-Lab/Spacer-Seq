cd /

for %%i in (01, 02, 03, 04, 05, 06, 07, 08, 09) do (
python "C:\Users\santi.bhattaraikline\Dropbox (Gladstone)\eVOLVER_and_Retro_Record/MiSeq_Data/SPCR_blast_v2.py" 2 %%i local BL21_DUETb3v2msdr_Cas_mphR86RT)

for /L %%i in (10,1,42) do (
python "C:\Users\santi.bhattaraikline\Dropbox (Gladstone)\eVOLVER_and_Retro_Record/MiSeq_Data/SPCR_blast_v2.py" 2 %%i local BL21_DUETb3v2msdr_Cas_mphR86RT)

cmd /k
