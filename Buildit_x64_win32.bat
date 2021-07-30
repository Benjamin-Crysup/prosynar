IF NOT EXIST builds MKDIR builds
IF EXIST builds\bin_x64_win32 RMDIR /S /Q builds\bin_x64_win32
IF EXIST builds\obj_x64_win32 RMDIR /S /Q builds\obj_x64_win32
IF EXIST builds\inc_x64_win32 RMDIR /S /Q builds\inc_x64_win32

MKDIR builds\bin_x64_win32
MKDIR builds\obj_x64_win32
MKDIR builds\inc_x64_win32

COPY includes_com\*.h builds\inc_x64_win32
COPY includes_prosynar\*.h builds\inc_x64_win32
FOR %%v in (source_com\*.cpp) DO g++ -mthreads %%v -g -O0 -Wall -Ibuilds\inc_x64_win32 -c -o builds\obj_x64_win32\%%~nv.o
FOR %%v in (source_com_any_win32\*.cpp) DO g++ -mthreads %%v -g -O0 -Wall -Ibuilds\inc_x64_win32 -c -o builds\obj_x64_win32\%%~nv.o
FOR %%v in (source_com_x64_any\*.cpp) DO g++ -mthreads %%v -g -O0 -Wall -Ibuilds\inc_x64_win32 -c -o builds\obj_x64_win32\%%~nv.o
FOR %%v in (source_prosynar\*.cpp) DO g++ -mthreads %%v -g -O0 -Wall -Ibuilds\inc_x64_win32 -c -o builds\obj_x64_win32\%%~nv.o
g++ -mthreads -o builds\bin_x64_win32\prosynar.exe builds\obj_x64_win32\*.o -lz -static-libgcc -static-libstdc++ -static -lpthread

builds\bin_x64_win32\prosynar --helpdumpgui | python utilities\ArgGuiBuilder.py > builds\bin_x64_win32\prosynargui.py
builds\bin_x64_win32\prosynar -- Fpprobreg --helpdumpgui | python utilities\ArgGuiBuilder.py > builds\bin_x64_win32\prosynarFpprobreggui.py
builds\bin_x64_win32\prosynar -- Fprover --helpdumpgui | python utilities\ArgGuiBuilder.py > builds\bin_x64_win32\prosynarFprovergui.py
builds\bin_x64_win32\prosynar -- Frover --helpdumpgui | python utilities\ArgGuiBuilder.py > builds\bin_x64_win32\prosynarFrovergui.py
builds\bin_x64_win32\prosynar -- Malign --helpdumpgui | python utilities\ArgGuiBuilder.py > builds\bin_x64_win32\prosynarMaligngui.py
builds\bin_x64_win32\prosynar -- Mflash --helpdumpgui | python utilities\ArgGuiBuilder.py > builds\bin_x64_win32\prosynarMflashgui.py
builds\bin_x64_win32\prosynar -- Mpear --helpdumpgui | python utilities\ArgGuiBuilder.py > builds\bin_x64_win32\prosynarMpeargui.py
COPY GUI.py builds\bin_x64_win32\

builds\bin_x64_win32\prosynar --helpdumpgui | python utilities\ManPageBuilder.py prosynar "merge paired end reads" prosynarmandesc.txt > builds\bin_x64_win32\prosynar.1
builds\bin_x64_win32\prosynar -- Fpprobreg --helpdumpgui | python utilities\ManPageBuilder.py prosynar_Fpprobreg "filter pairs with overlap in repeat" prosynarFpprobregmandesc.txt > builds\bin_x64_win32\prosynar_Fpprobreg.1
builds\bin_x64_win32\prosynar -- Fprover --helpdumpgui | python utilities\ManPageBuilder.py prosynar_Fprover "filter pairs with insufficient reference overlap" prosynarFprovermandesc.txt > builds\bin_x64_win32\prosynar_Fprover.1
builds\bin_x64_win32\prosynar -- Frover --helpdumpgui | python utilities\ManPageBuilder.py prosynar_Frover "filter pairs with insufficient naive reference overlap" prosynarFrovermandesc.txt > builds\bin_x64_win32\prosynar_Frover.1
builds\bin_x64_win32\prosynar -- Malign --helpdumpgui | python utilities\ManPageBuilder.py prosynar_Malign "merge pairs by alignment" prosynarMalignmandesc.txt > builds\bin_x64_win32\prosynar_Malign.1
builds\bin_x64_win32\prosynar -- Mflash --helpdumpgui | python utilities\ManPageBuilder.py prosynar_Mflash "merge pairs the way FLASH does it" prosynarMflashmandesc.txt > builds\bin_x64_win32\prosynar_Mflash.1
builds\bin_x64_win32\prosynar -- Mpear --helpdumpgui | python utilities\ManPageBuilder.py prosynar_Mpear "merge pairs the way PEAR does it" prosynarMpearmandesc.txt > builds\bin_x64_win32\prosynar_Mpear.1

