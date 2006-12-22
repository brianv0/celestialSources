@echo off
if .%1==. (set tag=Win32Debug ) else set tag=%1
set tempfile=%HOME%\tmpsetup.bat
d:\ground\CMT\v1r10p20011126\VisualC\cmt.exe -quiet -bat -pack=GRB -version=v0 setup -tag=%tag% >"%tempfile%"
if exist "%tempfile%" call "%tempfile%"
if exist "%tempfile%" del "%tempfile%"
set PATH=%LD_LIBRARY_PATH%;%PATH%