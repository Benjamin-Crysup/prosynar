mkdir -p builds
rm -rf builds/bin_x64_linux
rm -rf builds/obj_x64_linux
rm -rf builds/inc_x64_linux

mkdir builds/bin_x64_linux
mkdir builds/obj_x64_linux
mkdir builds/inc_x64_linux

cp includes_com/*.h builds/inc_x64_linux/
cp includes_prosynar/*.h builds/inc_x64_linux/
for v in source_com/*.cpp; do g++ -pthread $v -g -O0 -Wall -Ibuilds/inc_x64_linux -c -o builds/obj_x64_linux/$(basename $v .cpp).o; done;
for v in source_com_any_linux/*.cpp; do g++ -pthread $v -g -O0 -Wall -Ibuilds/inc_x64_linux -c -o builds/obj_x64_linux/$(basename $v .cpp).o; done;
for v in source_com_x64_any/*.cpp; do g++ -pthread $v -g -O0 -Wall -Ibuilds/inc_x64_linux -c -o builds/obj_x64_linux/$(basename $v .cpp).o; done;
for v in source_prosynar/*.cpp; do g++ -pthread $v -g -O0 -Wall -Ibuilds/inc_x64_linux -c -o builds/obj_x64_linux/$(basename $v .cpp).o; done;
g++ -pthread -o builds/bin_x64_linux/prosynar builds/obj_x64_linux/*.o -lz -ldl

./builds/bin_x64_linux/prosynar --helpdumpgui | python3 utilities/ArgGuiBuilder.py > builds/bin_x64_linux/prosynargui.py
./builds/bin_x64_linux/prosynar -- Fpprobreg --helpdumpgui | python3 utilities/ArgGuiBuilder.py > builds/bin_x64_linux/prosynarFpprobreggui.py
./builds/bin_x64_linux/prosynar -- Fprover --helpdumpgui | python3 utilities/ArgGuiBuilder.py > builds/bin_x64_linux/prosynarFprovergui.py
./builds/bin_x64_linux/prosynar -- Frover --helpdumpgui | python3 utilities/ArgGuiBuilder.py > builds/bin_x64_linux/prosynarFrovergui.py
./builds/bin_x64_linux/prosynar -- Malign --helpdumpgui | python3 utilities/ArgGuiBuilder.py > builds/bin_x64_linux/prosynarMaligngui.py
./builds/bin_x64_linux/prosynar -- Mflash --helpdumpgui | python3 utilities/ArgGuiBuilder.py > builds/bin_x64_linux/prosynarMflashgui.py
./builds/bin_x64_linux/prosynar -- Mpear --helpdumpgui | python3 utilities/ArgGuiBuilder.py > builds/bin_x64_linux/prosynarMpeargui.py
cp GUI.py builds/bin_x64_linux/

./builds/bin_x64_linux/prosynar --helpdumpgui | python3 utilities/ManPageBuilder.py prosynar "merge paired end reads" prosynarmandesc.txt > builds/bin_x64_linux/prosynar.1
./builds/bin_x64_linux/prosynar -- Fpprobreg --helpdumpgui | python3 utilities/ManPageBuilder.py prosynar_Fpprobreg "filter pairs with overlap in repeat" prosynarFpprobregmandesc.txt > builds/bin_x64_linux/prosynar_Fpprobreg.1
./builds/bin_x64_linux/prosynar -- Fprover --helpdumpgui | python3 utilities/ManPageBuilder.py prosynar_Fprover "filter pairs with insufficient reference overlap" prosynarFprovermandesc.txt > builds/bin_x64_linux/prosynar_Fprover.1
./builds/bin_x64_linux/prosynar -- Frover --helpdumpgui | python3 utilities/ManPageBuilder.py prosynar_Frover "filter pairs with insufficient naive reference overlap" prosynarFrovermandesc.txt > builds/bin_x64_linux/prosynar_Frover.1
./builds/bin_x64_linux/prosynar -- Malign --helpdumpgui | python3 utilities/ManPageBuilder.py prosynar_Malign "merge pairs by alignment" prosynarMalignmandesc.txt > builds/bin_x64_linux/prosynar_Malign.1
./builds/bin_x64_linux/prosynar -- Mflash --helpdumpgui | python3 utilities/ManPageBuilder.py prosynar_Mflash "merge pairs the way FLASH does it" prosynarMflashmandesc.txt > builds/bin_x64_linux/prosynar_Mflash.1
./builds/bin_x64_linux/prosynar -- Mpear --helpdumpgui | python3 utilities/ManPageBuilder.py prosynar_Mpear "merge pairs the way PEAR does it" prosynarMpearmandesc.txt > builds/bin_x64_linux/prosynar_Mpear.1


