#!/usr/bin/python
import sys
import getopt
import re
import string

p=re.compile('^<pattern id="12" $');

dane=list();
plikD=open("data.dat","r");
for linia in plikD:
  dane.append(int(linia[0:len(linia)-1]));
#for linia in dane:
  #print linia;

def sortujRosnaco(dane):
  j=1;
  while j<len(dane):
    key=dane[j];
    i=j-1;
    while i>-1 and dane[i]>key:
      dane[i+1]=dane[i];
      i=i-1;
    dane[i+1]=key;
    j=j+1;

sortujRosnaco(dane);
plikIn=open("ParticleData.xml","r");
linia="";
dane2=list();
for linia in plikIn:
 dane2.append(linia[0:len(linia)-1]);
plikD.close();
plikIn.close()
rozmiar=len(dane2);
i=0;
wsk=0;
pattern="^<particle id=\"";

fopen=open("ParticleDataBase.xml","w");
fopen.write("<?xml version=\"1.0\" encoding=\"UTF-8\"?>"+"\n\n");
fopen.write("<particleDataBase xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\">"+"\n");

while i<rozmiar and wsk<len(dane):
  p=re.compile(pattern+str(dane[wsk])+"\"");
  if p.match(dane2[i]):
    fopen.write(dane2[i]+"\n");
    while dane2[i]!="</particle>":
      i=i+1;
      fopen.write(dane2[i]+"\n");
    wsk=wsk+1;
    fopen.write("\n");
  i=i+1;
fopen.write("</particleDataBase>");
fopen.close();

