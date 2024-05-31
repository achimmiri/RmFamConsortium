#!/usr/bin/python3
import zipfile
from os import path
import sqlite3
import pandas
import sys
from rich import console
from rich import table
from rich import progress
from plot import main as sunburst

###fileToCheck = "/home/achimmir/temp/RmFamConsortium/R/mFam-data-master-0568feb3cb2847c0625ea5abd4e73cb0b2be9307.1.zip"
### This is the Zip file of the all the RMassBank files of the contributors
fileToCheck = sys.argv[1]

### Program shold be run as python3 create.suburst.py fileToCheck

def analyseFile(input, file):
    with input.open(file) as f:

        p = path.split(path.split(file)[0])

        data = {"Folder": p[1]}
        for line in f:
            decoded = line.decode('utf-8')
            decoded = decoded.strip(" \n\r\t")

            fields = decoded.split(":")
            if len(fields) == 2:
                if fields[0] == "AUTHORS":
                    data["Authors"] = fields[1]
                elif fields[0] == "CH$LINK" and "INCHIKEY" in fields[1]:
                    inchi = fields[1].strip().split(" ")
                    data["Inchikey"] = inchi[1]
                elif fields[0] == "CH$LINK" and "ChemOnt" in fields[1]:
                    ontology = fields[1].strip()[8:].strip()
                    data["Ontology"] = ontology
                elif fields[0] == "AC$MASS_SPECTROMETRY" and "ION_MODE" in fields[1]:
                    mode = fields[1].strip().split(" ")
                    data["Ionmode"] = mode[1]
    if "Authors" in data:
        return [data["Folder"], data["Authors"], data["Inchikey"], data["Ionmode"], data["Ontology"]]
    else:
        return None

def find(df, folder):
    for (_,row) in df.iterrows():
        if row["Folder"] == folder:
            return row["Value"]
    return 0

def main():
    input = zipfile.ZipFile(fileToCheck)
    inputfiles = input.namelist()
    data = []
    for file in progress.track(inputfiles, "Loading Files..."):
        tmp = analyseFile(input, file)
        if tmp is not None:
            data.append(tmp)
    df = pandas.DataFrame(data, columns=["Folder", "Authors", "Inchikey", "Ionmode", "Ontology"])
    conn = sqlite3.connect('file:cachedb?mode=memory&cache=shared')
    df.to_sql("data", conn)

    count = pandas.read_sql("select Folder, count(*) as Value from data group by Folder order by Folder", conn)
    pos = pandas.read_sql("select Folder,count(*) as Value from data where Ionmode='POSITIVE' group by Folder order by Folder",conn)
    neg = pandas.read_sql("select Folder,count(*) as Value from data where Ionmode='NEGATIVE' group by Folder order by Folder",conn)
    ucount = pandas.read_sql("select Folder, count(distinct(Inchikey)) as Value from data group by Folder order by Folder", conn)
    upos = pandas.read_sql("select Folder,count(distinct(Inchikey)) as Value from data where Ionmode='POSITIVE' group by Folder order by Folder",conn)
    uneg = pandas.read_sql("select Folder,count(distinct(Inchikey)) as Value from data where Ionmode='NEGATIVE' group by Folder order by Folder",conn)

    tcount = pandas.read_sql("select count(*) as Value from data", conn)
    tpos = pandas.read_sql("select count(*) as Value from data where Ionmode='POSITIVE'", conn)
    tneg = pandas.read_sql("select count(*) as Value from data where Ionmode='NEGATIVE'", conn)
    utcount = pandas.read_sql("select count(distinct(Inchikey)) as Value from data", conn)
    utpos = pandas.read_sql("select count(distinct(Inchikey)) as Value from data where Ionmode='POSITIVE'", conn)
    utneg = pandas.read_sql("select count(distinct(Inchikey)) as Value from data where Ionmode='NEGATIVE'", conn)

    data = pandas.read_sql("select a.Ontology, count(a.Inchikey) as count from (select distinct Inchikey, Ontology from data) as a group by a.Ontology", conn)
    
    data1=pandas.read_sql("select Ontology, count(Inchikey) as count from data group by Ontology",conn)

    data.to_csv("for_plot.csv", sep="\t", columns=['Ontology', 'count'], index=False)
    data1.to_csv("for_plot1.csv",sep="\t", columns=['Ontology', 'count'], index=False)


    data = []
    for (_,row) in count.iterrows():
        folder = row["Folder"]
        tmp = [folder, row["Value"], find(pos, folder), find(neg, folder), find(ucount, folder), find(upos, folder), find(uneg, folder)]
        data.append(tmp)

    tab = table.Table(title="Spectra Counts")
    tab.add_column("Contributor Dataset", no_wrap=True)
    tab.add_column("Number of spectra")
    tab.add_column("Number of Positive Spectra")
    tab.add_column("Number of Negative Spectra")
    tab.add_column("Number of Unique Compounds")
    tab.add_column("Number of Unique Compounds Positive")
    tab.add_column("Number of Unique Compounds Negative")

    for row in data:
        tab.add_row(str(row[0]),str(row[1]),str(row[2]),str(row[3]),str(row[4]),str(row[5]),str(row[6]))

    tab2 = table.Table(title="Totals")
    tab2.add_column("mFam_Consortium_library")
    tab2.add_column("POS")
    tab2.add_column("NEG")
    tab2.add_column("Total")

    tab2.add_row("Total no of spectra", str(tpos["Value"][0]), str(tneg["Value"][0]), str(tcount["Value"][0]))
    tab2.add_row("Unique Compounds", str(utpos["Value"][0]), str(utneg["Value"][0]), str(utcount["Value"][0]))


    con = console.Console()
    con.print(tab)
    con.print(tab2)

    sunburst("for_plot.csv","unique_plot.png")

    sunburst("for_plot1.csv","total_plot.png")


if __name__ == "__main__":
    main()
