# fxtPicoAna
Codes for FXT picoDst Analysis

# 1. Setup/ Sourcing
Login to BNL RACF computer and get the souce code:
```
git clone https://github.com/cdxing/fxtPicoAna.git

```

# 2. Environment of Working on RCF
```
stardev
kinit
aklog
cvs co StRoot/StPicoEvent
```

# 3. Compilation and Execution
Inside the /StRoot/StPicoEvent, do

```
make
```
To run the macro
```
root4star -b -q RunAnalyzer.C+
root4star -b -q PicoAnalyzer3.C+'("raw_PicoDst.root","output_name",2)'
```
where the "2" is the order of flow

# 4.Updating The Package/Setup
```
git pull
```

inside the /StRoot/StPicoEvent
```
make clean
make
```
