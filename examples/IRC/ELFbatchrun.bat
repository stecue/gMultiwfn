for /f %%i in ('dir c:\DA_IRCchk\*.fch /b') do (
Multiwfn c:\DA_IRCchk\%%i < ELFbatch.txt
rename ELF.cub %%~ni.cub
)