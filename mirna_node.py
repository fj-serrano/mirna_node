# -*- coding: utf-8 -*-


import pandas as pd
import matplotlib.pyplot as plt
import os
import sys
import argparse
import seaborn as sns


# FUNCIÓN QUE PERMITE AL USUARIO SELECCIONAR ARCHIVOS Y GUARDAR LAS RUTAS EN UNA LISTA
def select_archives():
    files = []
    # Argumentos de entrada para el script que utiliza para localizar los archivos de entrada, 
    # el nombre del estudio, método de normalización y localización de salida de salida
    parser = argparse.ArgumentParser(description = "Procesa un archivo de entrada")
    parser.add_argument("-i", "--input", required = True, help = "Ruta del archivo de entrada")
    parser.add_argument("-o", "--output", required = True, help = "Ruta del archivo de salida")
    parser.add_argument("-s", "--study", help = "Elige el nombre del estudio", default = "SRP")
    parser.add_argument("-d", "--db_relations", help = "Elige la ruta del archivo de relaciones entre bases de datos", default = "/shared/bak/TFG/serrano/scripts_outputs/db_relations/db_relations.tsv")
    parser.add_argument("-sp", "--specie", required = True, help = "Especie a la que pertenecen los datos", default = "hsa")
    args = parser.parse_args()
    file_path = args.input
    output_path = args.output
    study = args.study
    db_relations = args.db_relations
    specie = args.specie.lower()
    # Normaliza el path para evitar problemas con rutas relativas o absolutas
    file_path_clean = os.path.normpath(os.path.abspath(file_path))
    output_path = os.path.normpath(os.path.abspath(output_path))
    # Busca todas las carpetas cuyo nombre empiece por "SRP" y, dentro de cada una,
    # busca recursivamente una subcarpeta llamada igual que el argumento 'regression'
    # Si la encuentra, añade todas las .tsv que empiecen por DESeq2, edgeR o limma.
    metodos_prefijos = ("DESeq2", "edgeR", "limma")
    for root, dirs, files_in_dir in os.walk(file_path_clean):
        # Identifica directorios de estudio seleccionado
        if os.path.basename(root).startswith(study):
            # Recorre las subcarpetas del estudio en busca de la carpeta de regresion
            for subroot, subdirs, subfiles in os.walk(root):
                if os.path.basename(subroot).startswith("de_"):
                    # Añade solo los archivos .tsv cuyos nombres empiecen por alguno de los metodos y en los que se enfrenten 2 condiciones (vs)
                    for file in subfiles:
                        if file.endswith(".tsv") and file.startswith(metodos_prefijos) and ("vs" in file):
                            file_path_full = os.path.join(subroot, file)
                            files.append(file_path_full)
    return files, output_path, db_relations, specie

# FUNCIÓN QUE LEE LOS DISTINTOS ARCHIVOS SELECCIONADOS EN DATA FRAMES, EXTRAER Y OPERAR LA INFORMACIÓN IMPORTANTE EN ELLOS
def files_to_df(files):
    if not files:
        print("No se han seleccionado archivos.")
        sys.exit(1)
    else:
        # Lista con los métodos
        metodos = ("DESeq2", "edgeR", "limma")
        # Lista con las bases de datos
        dbs = ["miRBase", "mirgenedb", "mircarta"]
        # Lee cada archivo seleccionado, guardando los miRNAs significativos y sus valores de log2FoldChange y p-value ajustado
        # y guardando la base de datos de la que provienen y el método de normalización
        data = []
        for file in files:
            # Lee el archivo y guarda todos los miRNAs con padj < 0.05 y log2FoldChange > 1 o < -1
            df = pd.read_csv(file, sep = "\t")
            df_filtrado = df[(df["padj"] < 0.05) & ((df["log2FoldChange"] > 1) | (df["log2FoldChange"] < -1))]
            mirna = df_filtrado["gene"].tolist()
            log2FoldChange = df_filtrado["log2FoldChange"].tolist()
            padj = df_filtrado["padj"].tolist()

            # Extrae la base de datos y el método del nombre del archivo
            metodo = "UNDEFINED"
            base_de_datos = "UNDEFINED"
            regression = "UNDEFINED"
            # Busca el método, la base de datos y el nombre del estudio en la ruta del archivo
            ruta = file.split(os.sep)
            for part in ruta:
                if part.startswith("SRP"):
                    # Si empieza con SRP se asume que es el nombre del estudio, sigue la estructura "SRP_DB"
                    estudio_db = part.split("_")
                    estudio = estudio_db[0]
                    if estudio_db[1] in dbs:
                        base_de_datos = estudio_db[1]
                if part.startswith(metodos):
                    # Si empieza con método se asume que es el método, sigue la estructura "METODO_....tsv"
                    metodo_etc = part.split("_")
                    metodo = metodo_etc[0]
                if part.startswith("de_"):
                    # Si empieza con de_ se asume que es la carpeta de regresión, sigue la estructura "de_COND1_vs_COND2"
                    regression = part
            # Agrega la información a data como una lista de diccionarios
            for i in range(len(mirna)):
                data.append({"mirna": mirna[i], "log2FoldChange": log2FoldChange[i], "padj": padj[i], "estudio": estudio, "db": base_de_datos, "metodo": metodo, "assignment": regression})
        # Crea un DataFrame a partir de la lista de diccionarios
        df = pd.DataFrame(data)
        df_sorted = df.sort_values(by = ["estudio", "db", "metodo", "assignment"]).reset_index(drop = True)
    return df_sorted
       
# FUNCIÓN QUE LEE EL RESULTADO DE "db_relations.py" Y CONTRASTA SI LOS MIRNAS SIGNIFICATIVOS EXISTEN
# EN LAS OTRAS BASES DE DATOS, Y SI LO HACEN AÑADE COLUMNAS CON SUS NOMBRES ALTERNATIVOS Y LOS NODOS
def mirna_relation(df, db_relations, specie):
    # Guardamos en un DataFrame la tabla de relaciones entre bases de datos
    df_db_relations = pd.read_csv(db_relations, sep = "\t")
    # Creamos listas vacias para contar el número de bases de datos en las que aparece cada miRNA y para guardar los nombres alternativos
    db_counts = []
    db_alternative = []
    db_nodes = []

    # Para cada miRNA significativo, se mira su base de datos, en esa columna de df_db_relations se busca el miRNA
    # y si se encuentra, se añaden a df una columna con el número de bases de datos en las que aparece y en otra columna los nombres alternativos de ese miRNA en las otras bases de datos
    for index, row in df.iterrows():
        # Se obtiene la base de datos del miRNA
        db = row["db"] + (f"_{specie}")

        # Se inicializan el contador de bases de datos y la cadena de nombres alternativos
        count = int(0)
        alternative_mirnas = []
        alternative_nodes = []

        # Se busca un match del miRNA en la columna de su base de datos en df_db_relations
        row_db_relations = df_db_relations[df_db_relations[db] == (row["mirna"])]
        # Si row_db_relations está vacia se vuelve a probar, pero añadiendo un asterisco al final del nombre del mirna
        if row_db_relations.empty:
            row_db_relations = df_db_relations[df_db_relations[db] == (f"{row["mirna"]}*")]

        # Si no está vacío, se recorren las columnas de esa fila, y si el valor no es "-", se incrementa el contador y se añade el nombre alternativo a la cadena de nombres alternativos 
        if not row_db_relations.empty:
            for columna in row_db_relations.columns.tolist():
                if columna != db:
                   if f"{specie}" in columna:
                       mirnas = []
                       for i in range(len(row_db_relations)):
                           match_row = row_db_relations.iloc[i]
                           if match_row[columna] != "-":
                               count += 1
                               mirnas.append(match_row[columna])
                       if len(mirnas) == 0:
                           mirnas.append("-")
                       string_mirnas = ",".join(mirnas)
                       alternative_mirnas.append(f"{columna}:{string_mirnas}")
                   else:
                       nodes = []
                       for i in range(len(row_db_relations)):
                           match_row = row_db_relations.iloc[i]
                           if match_row[columna] != "-":
                               nodes.append(match_row[columna])
                       if len(nodes) == 0:
                           nodes.append("-")
                       string_nodes = ",".join(nodes)
                       alternative_nodes.append(f"{columna}:{string_nodes}")
                       
        db_counts.append(count)
        db_alternative.append(";".join(alternative_mirnas))
        db_nodes.append(";".join(alternative_nodes))

    # Añadimos todos los resultados de golpe a df
    df["db_count"] = db_counts
    df["db_alternativos"] = db_alternative
    df["nodes"] = db_nodes
    # Rellenamos huecos vacios con "-"
    df_final = df.fillna("-")
    df_final.to_csv("mirna_relation.tsv", sep = "\t", index = False)
    return df_final

# FUNCIÓN QUE TOMA EL DATA FRAME MIRNA_RELATION Y RESUME LA INFORMACIÓN DE CADA ESTUDIO POR BASE DE DATOS Y MÉTODO ESTADÍSTICO
def resumir_mirna_relation(df):
    # Para cada estudio se crea un diccionario que contiene cada columna y valor
    filas = []
    for estudio in df["estudio"].unique():
        for db in df["db"].unique():
            for metodo in df["metodo"].unique():
                for assignment in df["assignment"].unique():
                    # Recortamos el df para que contenga solo información relevante a un estudio concreto
                    db_recortado = df[(df["estudio"] == estudio) & (df["db"] == db) & (df["metodo"] == metodo) & (df["assignment"] == assignment)]
                    # Obtenemos el número de miRNAs significativos para todas esas condiciones
                    db_recortado_num_mirnas = len(db_recortado)
                    # Luego repetimos, pero buscanso solo los miRNAs que poseen nodos
                    db_recortado_nodo = db_recortado[~db_recortado["nodes"].str.contains("-")]
                    # Contamos también el número de miRNAs que contienen nodo
                    db_recortado_num_nodes = len(db_recortado_nodo)
                    # Y podemos calcular el % de miRNAs con nodo asegurándonos antes de cambiar el valor de db_recortado_num_mirnas a 1 en caso de que sea 0
                    if db_recortado_num_mirnas != 0:
                        porcentaje_mirnas_nodo = (db_recortado_num_nodes / db_recortado_num_mirnas) * 100
                    else:
                        porcentaje_mirnas_nodo = None

                    # Por último añadimos el diccionario completo a filas
                    filas.append({"estudio":estudio, "db": db, "metodo": metodo, "assignment": assignment, "n_mirnas_sig": db_recortado_num_mirnas, "n_mirnas_node": db_recortado_num_nodes, "mirnas_node_percent": porcentaje_mirnas_nodo})

    # Creamos el dataframe con toda la información y se guarda en un archivo
    df_resumen = pd.DataFrame(filas)
    df_resumen = df_resumen.sort_values(by = ["estudio", "db", "metodo", "assignment"]).reset_index(drop = True)
    df_resumen.to_csv("mirna_relation_resumen.tsv", sep = "\t", index = False)
    return df_resumen

# FUNCIÓN QUE REPRESENTA GRÁFICAMENTE LA INFORMACIÓN RESUMIDA PARA UN FÁCIL ANÁLISIS
def representacion_grafica(df):
    sns.set_theme(style = "ticks")
    # Se realizan 3 subplots, uno por db y compartiendo el eje Y
    fig, axex = plt.subplots(1, 3, figsize = (15, 5.5), sharey = True)
    databases = df["db"].unique()

    # Elegimos los colores para las cajas y los puntos
    box_colors = {"de_rcadj": "#9ecae1", "de_rcsa": "#fdae6b"}
    point_colors = {"de_rcadj": "#3182bd", "de_rcsa": "#e6550d"}

    for i, db in enumerate(databases):
        ax = axex[i]
        df_sub = df[df["db"] == db]
        # Pintamos un boxplot con seaborn, con el método en el eje X, el porcentaje de miRNAs con nodo en el eje Y
        sns.boxplot(data = df_sub, x = "metodo", y = "mirnas_node_percent", hue = "assignment", ax = ax, palette = box_colors, width = 0.6, dodge = True, fliersize = 0, boxprops = dict(alpha = 0.6))
        # Pintamos un stripplot con seaborn, con el método en el eje X, el porcentaje de miRNAs con nodo en el eje Y, y cada punto representando un estudio concreto
        sns.stripplot(data = df_sub, x = "metodo", y = "mirnas_node_percent", hue = "assignment", ax = ax, palette = point_colors, dodge = True, jitter = 0.15, size = 5, alpha = 0.8, marker = "o")

        # Cada panel se titula con el nombre de la base de datos, y se añaden etiquetas a los ejes
        ax.set_title(db, fontweight = "bold", pad = 12)
        ax.set_xlabel("Metodo de normalizacion", labelpad = 10)
        ax.set_xticklabels(["DESeq2", "edgeR", "limma"])
        if i == 0:
            ax.set_ylabel("Porcentaje de miRNAs con nodo", labelpad = 10)
        else:
            # Así evitamos que se repitan las etiquetas en los paneles siguientes
            ax.set_ylabel("")
        # Se eliminan las leyendas de cada panel para evitar que se repitan
        ax.get_legend().remove()
    
    # Se añade una leyenda común
    handels, labels = axex[0].get_legend_handles_labels()

    # Tomamos solo los primeros 2 elementos (cajas) para la leyenda
    fig.legend(handels[:2], ["de_rcadj", "de_rcsa"], loc = "upper left", ncol = 2, title = "Asignacion de lecturas", title_fontsize = "11", fontsize = "9", bbox_to_anchor = (1.02, 1), frameon = True)
    sns.despine()
    plt.tight_layout()
    plt.savefig("importancia_biologica_mirnas.png", dpi = 300, bbox_inches = "tight")

    return

# FUNCIÓN QUE CREA UN DATAFRAME CON LOS MIRNAS DIFERENCIALMENTE EXPRESADOS EN MIRGENEDB QUE NO SE DETECTAN COMO DIFERENCIALMENTE EXPRESADOS EN MIRBASE Y MIRCARTA
def mirna_de_mirgenedb(df):
    data = []
    # Se copia el dataframe para no modificar el original
    df_nodes = df.copy()
    # Se recorta el dataframe para quedarnos solo con los miRNAs de MirGeneDB
    df_mirgenedb = df_nodes[df_nodes["db"] == "mirgenedb"]

    # Se recorre cada miRNA de df_mirgenedb, se toman sus nombres alternativos de la columna "db_alternativos" y se comprueba si aparecen como diferencialmente expresados 
    # en miRBase o MirCarta, y si no es así, se añaden a una lista de miRNAs diferencialmente expresados en MirGeneDB pero no en las otras bases de datos
    for index, row in df_mirgenedb.iterrows():
        # Primero se obtiene diccionario con los nombres alternativos del miRNA
        alternative_names = {}
        texto_separado = row["db_alternativos"].split(";")
        for item in texto_separado:
            db_and_name = item.split(":")
            alternative_names[db_and_name[0].split("_")[0]] = db_and_name[1]
            
        # Se comprueba si el miRNA aparece como diferencialmente expresado en miRBase o MirCarta
        for db in alternative_names.keys():
            alternative_mirnas = alternative_names[db].split(",")
            for alt_mirna in alternative_mirnas:
                if alt_mirna != "-":
                    df_alt = df_nodes[(df_nodes["mirna"] == alt_mirna) & (df_nodes["db"] == db) & (df_nodes["estudio"] == row["estudio"]) & (df_nodes["metodo"] == row["metodo"]) & (df_nodes["assignment"] == row["assignment"])]
                    if df_alt.empty:
                        data.append({"mirna_mirgenedb": row["mirna"], "mirna_alternativo": alt_mirna, "db_alternativo": db, "estudio": row["estudio"], "metodo": row["metodo"], "assignment": row["assignment"]})
                else:
                    data.append({"mirna_mirgenedb": row["mirna"], "mirna_alternativo": alt_mirna, "db_alternativo": db, "estudio": row["estudio"], "metodo": row["metodo"], "assignment": row["assignment"]})
    
    # Una vez se han recorrido todos los miRNAs, se crea un dataframe con la información recopilada y se guarda en un archivo
    df_mirna_de_mirgenedb = pd.DataFrame(data)
    df_mirna_de_mirgenedb = df_mirna_de_mirgenedb.sort_values(by = ["estudio", "metodo", "assignment", "db_alternativo"]).reset_index(drop = True)
    df_mirna_de_mirgenedb.to_csv("mirna_de_mirgenedb.tsv", sep = "\t", index = False)
    return df_mirna_de_mirgenedb

# FUNCIÓN QUE TOMA EL DATAFRAME MIRNA_DE_MIRGENEDB Y RESUME LA INFORMACIÓN
def resumir_mirna_de_mirgenedb(df):
    # Mezclamos los resultados de todos los estudios buscando los miRNAs que más se repiten respetando las diferencias por bases de datos
    df_mirna_de_mirgenedb = df.copy()
    # Se agrupa por el nombre del miRNA de MirGeneDB y se cuenta el número de veces que aparece cada miRNA
    df_resumen_mirna_de_mirgenedb = df_mirna_de_mirgenedb.groupby(["mirna_mirgenedb", "mirna_alternativo", "db_alternativo"]).size().reset_index(name = "count")
    # Se ordena el dataframe por la base de datos y el número de veces que aparece cada miRNA
    df_resumen_mirna_de_mirgenedb = df_resumen_mirna_de_mirgenedb.sort_values(by = ["db_alternativo", "count"], ascending = False).reset_index(drop = True)
    # Se guarda el dataframe en un archivo
    df_resumen_mirna_de_mirgenedb.to_csv("mirna_de_mirgenedb_resumen.tsv", sep = "\t", index = False)
    return

def main():
    files, output_path, db_relations, specie = select_archives()
    output_path_new = output_path
    os.makedirs(output_path_new, exist_ok = True)
    os.chdir(output_path_new)
    df = files_to_df(files)
    df_relations = mirna_relation(df, db_relations, specie)
    df_resumen = resumir_mirna_relation(df_relations)
    representacion_grafica(df_resumen)
    df_mirna_de_mirgenedb = mirna_de_mirgenedb(df_relations)
    resumir_mirna_de_mirgenedb(df_mirna_de_mirgenedb)

if __name__ == "__main__":
    main()