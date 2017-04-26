@echo off
setlocal enabledelayedexpansion
for /f  %%i in ('dir *.gjf /b') do (
echo Calculating %%i...
g09 %%~ni %%~ni&if !errorlevel!==0 (
echo Calculation of %%i finished at !date! !time!) else (
echo Calculation of %%i failed at !date! !time!))