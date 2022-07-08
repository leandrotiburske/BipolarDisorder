`.raw` mass espectometry files were converted to `.mzXML` fomat using the following command:

```
for file in *.raw;
do docker run -it --rm -e WINEDEBUG=-all -v ~/Documents/IC/Lipidome/Data:/data chambm/pwiz-skyline-i-agree-to-the-vendor-licenses wine msconvert $file --mzXML;
done

```
