// Copyright [2018] [Pengyu Hong at Brandeis University]
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//     http://www.apache.org/licenses/LICENSE-2.0
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

public enum CMonosaccharideSet {
    // set of Monosaccharide
    Xyl("Xyl", "Xylose", "Xyl", 1,
            new int[]{3, 4, 5},
            new int[]{2, 3, 4, 5, 6},
            new CMonosaccharideInfo(164.0684734957,"C6-H12-O5",
                    new String[]{"O", "CH2O", "CH2O", "CH2O", "CH2O", "CH2"},
                    new double[]{15.9949146221, 30.0105646863, 30.0105646863, 30.0105646863, 30.0105646863, 14.0156500642}),
            new CMonosaccharideInfo(206.1154236883, "C8-H16-O5",
                    new String[]{"O", "C2H4O", "C2H4O", "C2H4O", "C2H4O", "CH2" },
                    new int[]{0, 0, 1, 1, 1, 0},
                    new double[]{15.9949146221, 44.0262147505, 44.0262147505, 44.0262147505, 44.0262147505, 14.0156500642})
    ),
    dHex(  "dHex","Fucose","Fuc", 2,
            new int[] {3, 4, 5},
            new int[] {2, 3, 4, 5, 6, 6},
            new CMonosaccharideInfo( 164.0684734957, "C6-H12-O5",
                    new String[]{"O", "CH2O", "CH2O", "CH2O", "CH2O", "C2H4"} ,
                    new double[]{15.9949146221,30.0105646863,30.0105646863,30.0105646863,30.0105646863,28.0313001284} ),
            new CMonosaccharideInfo( 220.1310737525, "C9-H18-O5",
                    new String[]{"O", "C2H4O", "C2H4O", "C2H4O", "C2H4O", "C2H4"} ,
                    new int[]{0, 0, 1, 1, 1, 0},
                    new double[]{15.9949146221,44.0262147505,44.0262147505,44.0262147505,44.0262147505,28.0313001284} ) ),
    Glc(  "Glc","Glucose","Hex", 3,
            new int[] {3, 4, 5, 6},
            new int[] {2, 3, 4, 5, 6, 6},
            new CMonosaccharideInfo( 180.0633881178, "C6-H12-O6",
                    new String[]{"O", "CH2O", "CH2O", "CH2O", "CH2O", "C2H4O"} ,
                    new double[]{15.9949146221,30.0105646863,30.0105646863,30.0105646863,30.0105646863,44.0262147505} ),
            new CMonosaccharideInfo( 250.1416384388, "C10-H20-O6",
                    new String[]{"O", "C2H4O", "C2H4O", "C2H4O", "C2H4O", "C3H6O"} ,
                    new int[]{0, 0, 1, 1, 1, 1},
                    new double[]{15.9949146221,44.0262147505,44.0262147505,44.0262147505,44.0262147505,58.0418648147} ) ),
    Gal(  "Gal","Galactose","Hex", 3,
            new int[] {3, 4, 5, 6},
            new int[] {2, 3, 4, 5, 6, 6},
            new CMonosaccharideInfo( 180.0633881178, "C6-H12-O6",
                    new String[]{"O", "CH2O", "CH2O", "CH2O", "CH2O", "C2H4O"} ,
                    new double[]{15.9949146221,30.0105646863,30.0105646863,30.0105646863,30.0105646863,44.0262147505} ),
            new CMonosaccharideInfo( 250.1416384388, "C10-H20-O6",
                    new String[]{"O", "C2H4O", "C2H4O", "C2H4O", "C2H4O", "C3H6O"} ,
                    new int[]{0, 0, 1, 1, 1, 1},
                    new double[]{15.9949146221,44.0262147505,44.0262147505,44.0262147505,44.0262147505,58.0418648147} ) ),
    Man(  "Man","Mannose","Hex", 3,
            new int[] {3, 4, 5, 6},
            new int[] {2, 3, 4, 5, 6, 6},
            new CMonosaccharideInfo( 180.0633881178, "C6-H12-O6",
                    new String[]{"O", "CH2O", "CH2O", "CH2O", "CH2O", "C2H4O"} ,
                    new double[]{15.9949146221,30.0105646863,30.0105646863,30.0105646863,30.0105646863,44.0262147505} ),
            new CMonosaccharideInfo( 250.1416384388, "C10-H20-O6",
                    new String[]{"O", "C2H4O", "C2H4O", "C2H4O", "C2H4O", "C3H6O"} ,
                    new int[]{0, 0, 1, 1, 1, 1},
                    new double[]{15.9949146221,44.0262147505,44.0262147505,44.0262147505,44.0262147505,58.0418648147} ) ),
    GlcA(  "GlcA","Glucuronic acid","HexA", 4,
            new int[] {3, 4, 5},
            new int[] {2, 3, 4, 5, 6, 6},
            new CMonosaccharideInfo( 194.0426526757, "C6-H10-O7",
                    new String[]{"O", "CH2O", "CH2O", "CH2O", "CH2O", "C2H2O2"} ,
                    new double[]{15.9949146221,30.0105646863,30.0105646863,30.0105646863,30.0105646863,58.0054793084} ),
            new CMonosaccharideInfo( 264.1209029967, "C10-H18-O7",
                    new String[]{"O", "C2H4O", "C2H4O", "C2H4O", "C2H4O", "C3H4O2"} ,
                    new int[]{0, 0, 1, 1, 1, 1},
                    new double[]{15.9949146221,44.0262147505,44.0262147505,44.0262147505,44.0262147505,72.0211293726} ) ),
    IdoA(  "IdoA","Iduronic acid","HexA", 4,
            new int[] {3, 4, 5},
            new int[] {2, 3, 4, 5, 6, 6},
            new CMonosaccharideInfo( 194.0426526757, "C6-H10-O7",
                    new String[]{"O", "CH2O", "CH2O", "CH2O", "CH2O", "C2H2O2"} ,
                    new double[]{15.9949146221,30.0105646863,30.0105646863,30.0105646863,30.0105646863,58.0054793084} ),
            new CMonosaccharideInfo( 264.1209029967, "C10-H18-O7",
                    new String[]{"O", "C2H4O", "C2H4O", "C2H4O", "C2H4O", "C3H4O2"} ,
                    new int[]{0, 0, 1, 1, 1, 1},
                    new double[]{15.9949146221,44.0262147505,44.0262147505,44.0262147505,44.0262147505,72.0211293726} ) ),
    GalNAc(  "GalNAc","N-acetylgalactosamine","HexNAc", 5,
            new int[] {4, 5, 6},
            new int[] {2, 3, 4, 5, 6, 6},
            new CMonosaccharideInfo( 221.0899372193, "C8-H15-N-O6",
                    new String[]{"O", "CH2O", "C3H5NO", "CH2O", "CH2O" ,"C2H4O"} ,
                    new double[]{15.9949146221,30.0105646863,71.0371137878,30.0105646863,30.0105646863,44.0262147505} ),
            new CMonosaccharideInfo( 291.1681875403, "C12-H23-N-O6",
                    new String[]{"O", "C2H4O", "C4H7NO", "C2H4O", "C2H4O", "C3H6O"} ,
                    new int[]{0, 0, 1, 1, 1, 1},
                    new double[]{15.9949146221,44.0262147505,85.0527638520,44.0262147505,44.0262147505,58.0418648147} ) ),
    GlcNAc(  "GlcNAc","N-acetylglucosamine","HexNAc", 5,
            new int[] {4, 5, 6},
            new int[] {2, 3, 4, 5, 6, 6},
            new CMonosaccharideInfo( 221.0899372193, "C8-H15-N-O6",
                    new String[]{"O", "CH2O", "C3H5NO", "CH2O", "CH2O" ,"C2H4O"} ,
                    new double[]{15.9949146221,30.0105646863,71.0371137878,30.0105646863,30.0105646863,44.0262147505} ),
            new CMonosaccharideInfo( 291.1681875403, "C12-H23-N-O6",
                    new String[]{"O", "C2H4O", "C4H7NO", "C2H4O", "C2H4O", "C3H6O"} ,
                    new int[]{0, 0, 1, 1, 1, 1},
                    new double[]{15.9949146221,44.0262147505,85.0527638520,44.0262147505,44.0262147505,58.0418648147} ) ),
    Kdo(  "Kdo","Kdo","Kdo", 6,
            new int[] {4, 5, 6},
            new int[] {2, 2, 3, 4, 5, 6, 6, 6},
            new CMonosaccharideInfo( 238.0688674262, "C8-H14-O8",
                    new String[]{"O", "C2H2O3", "CH2", "CH2O", "CH2O", "C3H6O2"} ,
                    new double[]{15.9949146221,74.0003939305,14.0156500642,30.0105646863,30.0105646863,74.0367794368} ),
            new CMonosaccharideInfo( 322.1627678114, "C14-H26-O8",
                    new String[]{"O", "C4H6O3", "CH2", "C2H4O", "C2H4O", "C5H10O2"} ,
                    new int[]{0, 1, 0, 1, 1, 2},
                    new double[]{15.9949146221,102.0316940589,14.0156500642,44.0262147505,44.0262147505,102.0680795652} ) ),
    NeuAc(  "NeuAc","N-acetylneuraminic acid","NeuAc", 7,
            new int[] {4, 6},
            new int[] {2, 2, 3, 4, 5, 6, 6, 6, 6},
            new CMonosaccharideInfo( 309.1059812140, "C11-H19-N-O9",
                    new String[]{"O", "C2H2O3", "CH2", "CH2O", "C3H5NO", "C4H8O3"} ,
                    new double[]{15.9949146221,74.0003939305,14.0156500642,30.0105646863,71.0371137878,104.0473441231} ),
            new CMonosaccharideInfo( 407.2155316634, "C17-H31-N-O9",
                    new String[]{"O", "C4H6O3", "CH2", "C2H4O", "C4H7NO", "C7H14O3"} ,
                    new int[]{0, 1, 0, 1, 1, 3},
                    new double[]{15.9949146221,102.0316940589,14.0156500642,44.0262147505,85.0527638520,146.0942943157} ) ),
    NeuGc(  "NeuGc","N-glycolylneuraminic acid","NeuGc", 8,
            new int[] {4, 6},
            new int[] {2, 2, 3, 4, 5, 6, 6, 6, 6},
            new CMonosaccharideInfo( 325.1008958361, "C11-H19-N-O10",
                    new String[]{"O", "C2H2O3", "CH2", "CH2O", "C3H5NO2", "C4H8O3"} ,
                    new double[]{15.9949146221,74.0003939305,14.0156500642,30.0105646863,87.0320284099,104.0473441231} ),
            new CMonosaccharideInfo( 437.2260963497, "C17-H31-N-O10",
                    new String[]{"O", "C4H6O3", "CH2", "C2H4O", "C5H9NO2", "C7H14O3"} ,
                    new int[]{0, 1, 0, 1, 1, 3},
                    new double[]{15.9949146221,102.0316940589,14.0156500642,44.0262147505,115.0633285383,146.0942943157} ) );


    final String symbol;
    final String name;
    final String sacClass;
    final int classID;
    final int[] legalLinkedInVs;
    final int[] carbon2vertex;
    final CMonosaccharideInfo sacNative;
    final CMonosaccharideInfo sacPermethylated;

    CMonosaccharideSet(String symbol, String name, String sacClass, int classID, int[] legalLinkedInVs, int[] carbon2vertex, CMonosaccharideInfo sacNative, CMonosaccharideInfo sacPermethylated) {
        this.symbol = symbol;
        this.name = name;
        this.sacClass = sacClass;
        this.classID = classID;
        this.legalLinkedInVs = legalLinkedInVs;
        this.carbon2vertex = carbon2vertex;
        this.sacNative = sacNative;
        this.sacPermethylated = sacPermethylated;
    }

    // inner static class for abstracting information of both permethylated and non-permed monosacs, must be static otherwise ENUM type values cannot initiate them inside ().
    static class CMonosaccharideInfo{
        final double mass;
        final String compostition;
        final String[] vertexComposition;
        final int[] vertexPermethylated;
        final double[] vertexMass;

        private CMonosaccharideInfo(double mass, String compostition, String[] vertexComposition, int[] vertexPermethylated, double[] vertexMass) {
            this.mass = mass;
            this.compostition = compostition;
            this.vertexComposition = vertexComposition;
            this.vertexPermethylated = vertexPermethylated;
            this.vertexMass = vertexMass;
        }

        private CMonosaccharideInfo(double mass, String compostition, String[] vertexComposition, double[] vertexMass) {
            this(mass, compostition, vertexComposition, null, vertexMass); // seems must use null to initiate, otherwise error
        }
    }
}
