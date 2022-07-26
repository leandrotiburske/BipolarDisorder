# Lipidome

## Raw files conversion

`.raw` mass espectometry files were converted to `.mzXML` fomat using the following Linux command:

```
for file in *.raw;
do docker run -it --rm -e WINEDEBUG=-all -v ~/yourfilespath/Data:/data chambm/pwiz-skyline-i-agree-to-the-vendor-licenses wine msconvert $file --mzXML;
done

```
