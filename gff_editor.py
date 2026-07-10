import re
from collections import defaultdict

# 1) CDS에서 gene별 Note 수집
notes = {}
matches_protein_ref = {}
valid_orf = {}
missing_stop_codon = {}


with open('Genome/Po1f.gff3') as f:
    for line in f:
        if line.startswith('#'): continue
        c = line.rstrip('\n').split('\t')
        
        if len(c) < 9: continue
        if c[2] == 'CDS':

            lt = re.search(r'locus_tag=([^;]+)', c[8])
            nt = re.search(r'Note=([^;]+)', c[8])

            if lt and nt and lt.group(1) not in notes:
                notes[lt.group(1)] = nt.group(1)
        
        elif c[2] == "mRNA":
            lt = re.search(r'locus_tag=([^;]+)', c[8])
            mpr = re.search(r'matches_protein_ref=([^;]+)', c[8])
            voc = re.search(r'valid_orf=([^;]+)', c[8])
            msc = re.search(r'missing_stop_codon=([^;]+)', c[8])

            if lt:
                matches_protein_ref[lt.group(1)] = mpr.group(1) if mpr else "False"
                valid_orf[lt.group(1)] = voc.group(1) if voc else "False"
                missing_stop_codon[lt.group(1)] = msc.group(1) if msc else "False"


# 2) gene 줄에 Note + description 주입
out = open('Genome/Po1f_annot.gff3', 'w')
with open('Genome/Po1f.gff3') as f:
    for line in f:
        if line.startswith('#'):
            out.write(line); continue
        c = line.rstrip('\n').split('\t')
        if len(c) >= 9 and c[2] == 'gene':

            lt = re.search(r'locus_tag=([^;]+)', c[8])
            if lt and lt.group(1) in notes:
                note = notes[lt.group(1)].replace('%2C', ',')
                # Name을 사람이 읽을 수 있게 교체
                

                status_text = ""
                if valid_orf.get(lt.group(1)) == "False":
                    status_text += "invalid_orf"
                elif missing_stop_codon.get(lt.group(1)) == "True":
                    status_text += "no_stop_codon"
                elif matches_protein_ref.get(lt.group(1)) == "False":
                    status_text += "no_protein_ref"

                if len(status_text) > 0:
                    status_text = f'({status_text}) '


                type_text = ""
                if note.find("weakly similar to uniprot|") != -1:
                    short = note.split("weakly similar to uniprot|")[1]
                    type_text = "uni_weak"

                elif note.find("similar to uniprot|") != -1:
                    short = note.split("similar to uniprot|")[1]
                    type_text = "uni"

                elif note.find("weakly similar to wi|") != -1:
                    short = note.split("weakly similar to wi|")[1]
                    type_text = "wi_weak"

                elif note.find("similar to wi|") != -1:
                    short = note.split("similar to wi|")[1]
                    type_text = "wi"

                elif note.find("similar ") != -1:
                    short = note.split("similar ")[1]
                    type_text = "similar"

                elif note.find("no similarity") != -1:
                    short = "K no similarity"

                else:
                    short = "K something went wrong!"
                
                short.replace("%2C", ",").replace("%20", " ").replace("%28", "(").replace("%29", ")").replace("%2F", "/").replace("%3A", ":").replace("%3B", ";").replace("%3D", "=").replace("%3F", "?").replace("%5B", "[").replace("%5D", "]").replace("%5E", "^").replace("%7B", "{").replace("%7C", "|").replace("%7D", "}").replace("%7E", "~")

                short_cut = re.search(r'([^;%,]+)', short).group(1)

                short_cut = short_cut[short_cut.find(" ")+1:]
                if len(short_cut) == 0:
                    print(f"Warning: {lt.group(1)} has no short cut!")
                    print(f"Note: {note}")
                    print(f"Short: {short}")

                if short_cut[0].capitalize() == short_cut[0]:
                    short_with_species = short_cut[0] + ". " + short_cut[short_cut.find(" ")+1:] 

                print(short_with_species)


                c[8] += f';Note={status_text}{note};description={status_text}{short_with_species}{type_text}'
                line = '\t'.join(c) + '\n'
        out.write(line)
out.close()



"""
    쓰이는 요소: Name, description 두 가지밖에 없음...

"""