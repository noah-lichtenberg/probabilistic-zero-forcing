import pandas as pd
import networkx as nx

# 1. Load Excel file. Gathered from the BEA website.
df = pd.read_excel("BEAData.xlsx")

# 2. Extract sector names (row 5, columns C..Q)
sector_names = df.iloc[5, 2:17].tolist()   # length 15

# 3. Build a mapping from int → sector name
id_to_sector = {i: name for i, name in enumerate(sector_names)}

# 4. Extract the 15×15 matrix (rows 8..22, cols C..Q)
# Excel rows 8–22  → pandas rows 7:22
# Excel columns C–Q → pandas cols 2:17
M = df.iloc[6:21, 2:17].replace('---', 0).astype(float).values

print(M)
# 5. Build directed graph with integer nodes
BEAinputoutputdata = nx.DiGraph()
BEAinputoutputdata.add_nodes_from(range(len(sector_names)))  # nodes = 0..14

# 6. Add weighted edges: (col j → row i), INCLUDING self-edges
n = len(sector_names)
for i in range(n):          # buyer (row)
    for j in range(n):      # seller (column)
        w = M[i, j]
        if w > 0:
            BEAinputoutputdata.add_edge(j, i, weight=w)

print("Self-loops:", [(u, v) for u, v in BEAinputoutputdata.edges() if u == v])