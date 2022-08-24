Etapes pour le debug :

1)

python3 filter.py <input_path> <output_path> 
    -> transforme un bucket kff validate en liste de Kmers type KMC_dump

2)

bash debug.sh <input_path>
    -> lance kff instr, bucket et compact -s avec k = 25 et m = 10
    + validate les outputs avec verbose dans fichier .txt

3) (option)

python3 verif.py <input_path>
    -> affiche le nombre de Kmer dans un fichier .txt obtenu par étape 2.
