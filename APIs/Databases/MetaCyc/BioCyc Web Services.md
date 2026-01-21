BioCyc Web Services
BioCyc offers several classes of web services that are summarized in the table below. These services include:
Data retrieval services that return data in XML, JSON, and/or column-delimited formats (BioCyc data can also be downloaded in a variety of XML and non-XML formats)

Visualization services, e.g., to map gene expression data onto pathway diagrams

Services for manipulating SmartTables
Use of the BioCyc web services API is subject to the terms of the BioCyc Databases Limited Use License. Note that to avoid congestion on the BioCyc web servers, we respectfully request that when writing programs that issue large numbers of web service queries, you limit requests to on average not more than one per second.
Note that in order to access BioCyc web services you must establish a session (see below) and utilize the returned cookies in subsequent web-services calls.


Class of Service	Service	Input	Output	HTTP GET or POST
Retrieve a Single PGDB Object	BioCyc object retrieval
BioCyc Object ID	ptools-XML	GET
BioPAX Pathway Data	BioCyc Pathway ID	BioPAX XML	GET
Search and Retrieve a Set of Objects
Retrieve objects by name or accession	BioCyc Object Name or Accession	tab delimited, JSON, ptools-XML	GET, POST
Foreign object-ID (external database identifier) retrieval	Foreign database ID	tab delimited, JSON, ptools-XML	GET, POST
Retrieve object sets, e.g., genes-of-pathway, compounds-of-pathway	BioCyc Object ID	ptools-XML	GET
Evaluate query written in BioVelo query language	BioVelo Query	Object set encoded as ptools-XML	
Metabolite monoisotopic weight retrieval	Monoisotopic MW, tolerance	tab delimited, JSON, ptools-XML	GET, POST
Metabolite chemical formula retrieval	Chemical Formula	tab delimited, JSON, ptools-XML	GET, POST
Metabolite Translation Service		tab delimited	GET, POST
Retrieve metabolites by SMILES string	SMILES String	tab delimited, JSON	GET, POST
Retrieve metabolites by INCHI string	INCHI String	tab delimited, JSON, ptools-XML	GET, POST
Visualization Services	Generate Compound Image	Compound ID	GIF image	GET
Generate Pathway Diagram	Pathway ID	GIF image	GET, POST
Genome Browser Visualization	Organism and Replicon Names and / or IDs	HTML	
Omics Visualization Services	Pathway Diagram Painted with Omics Data	Pathway ID	HTML	GET, POST
Table of Pathway Diagrams Painted with Omics Data	Pathway ID	HTML	GET, POST
Highlight/Paint Objects on Metabolic Network Diagram	Object Names and / or IDs	HTML	GET
Paint Single Omics Data from a File on Metabolic Network Diagram	Object Names and / or IDs	HTML	POST
Paint Multi Omics Data From File via POST	Multi Omics Text File	URL	POST
SmartTables
SmartTable Creation			PUT
SmartTable Data Retrieval		tab delimited, JSON, tsv	GET
SmartTable Add Transform Column			GET
SmartTable Add Property Column			GET
SmartTable Delete Row			GET
SmartTable Delete Column			GET
SmartTable Modify Name			GET
SmartTable Modify Description			GET
SmartTable Copy			GET
SmartTable Delete			GET
SmartTable Enrichment Analysis			GET
Establishing a Web Services session
To use the BioCyc web services, you must establish a session and log in to your BioCyc account by POSTing your email and password to https://websvc.biocyc.org/credentials/login/. The cookies returned from this request must then be passed back to BioCyc for every subsequent web service request.


Example using curl
    # Log in and save cookies to file biocyc_cookies (POST only):
    curl -c biocyc_cookies -d "email=[email]&password=[password]" -X POST https://websvc.biocyc.org/credentials/login/
    # Issue web service request, passing cookie file:
    curl -b biocyc_cookies https://websvc.biocyc.org/[request]
  
Example using python Requests module
    s = requests.Session() # create session
    # Post login credentials to session:
    s.post('https://websvc.biocyc.org/credentials/login/', data={'email':'[email]', 'password':'[password]'})
    # Issue web service request:
    r = s.get('https://websvc.biocyc.org/[request]')
  
Data Retrieval Web Services: Pathway Tools XML
Any pathway, reaction, compound, gene, protein, RNA or transcription-unit object in a BioCyc database can be retrieved in ptools-xml format, an XML format that is based on and closely resembles the underlying Pathway Tools schema. A single object can be requested using its BioCyc identifier, or a query can be issued using either a subset of the Pathway Tools API functions or the BioVelo query language to retrieve multiple objects.
The following documents will be useful in interpreting the ptools-xml format:

The Pathway Tools Schema Guide	A guide to the internal Pathway Tools representation, describing classes and their slots (attributes).
Guide to ptools-xml	This document describes the ptools-xml format and the differences between ptools-xml and the Pathway Tools Schema, including the mapping of class and slot names from one to the other.
ptools-xml.xsd	An XMLSchema document describing the ptools-xml format
ptools.wsdl	A WSDL document describing the ptools-xml-based web services.
Retrieve a Single BioCyc Object
Biocyc KB-Version Retrieval
This service finds the version of the KB and returns the result in JSON format.
The URLs for retrieving the KB-version is

https://websvc.biocyc.org/kb-version?orgid=[ORGID]

where

[ORGID] is the identifier for the organism database, e.g. ECOLI, META, BSUB.
Example URLs:

https://websvc.biocyc.org/kb-version?orgid=ECOLI
The default output is a JSON format as follows:

{"kb-version":"19.5"}


Biocyc Object-Id-Based Retrieval of XML Data
The URL to retrieve a single BioCyc object in ptools-xml format is

https://websvc.biocyc.org/getxml?[ORGID]:[OBJECT-ID] or
https://websvc.biocyc.org/getxml?id=[ORGID]:[OBJECT-id]&detail=[none|low|full]

where

[ORGID] is the identifier for the organism database, e.g. ECOLI, META, AFER243159.
[OBJECT-ID] is in the BioCyc identifier for the object, e.g. ARGSYN-PWY, EG11025, FRUCTOSE-6P. Note that object identifiers are case-sensitive.
[none|low|full] indicates whether the returned output should contain no detail, low detail or full detail for the requested object. If no detail parameter is supplied, the request defaults to full detail.
Example URLs:

https://websvc.biocyc.org/getxml?BSUB:GLYCOLYSIS
Retrieve the glycolysis pathway in Bacillus subtilis in ptools-xml format.
https://websvc.biocyc.org/getxml?id=META:ASPARTATEKIN-RXN
Retrieve the aspartate kinase reaction from MetaCyc in ptools-xml format.
https://websvc.biocyc.org/getxml?id=ECOLI:TRYPSYN-APROTEIN&detail=low
Retrieve the trpA gene product from EcoCyc at low detail level.
BioPAX Pathway Data
Pathway data for an individual pathway is available in BioPAX XML format (both BioPAX Level 2 and Level 3). The URL to access a pathway in BioPAX format is:
https://websvc.biocyc.org/[ORGID]/pathway-biopax?type=[2|3]&object=[PATHWAY]

where

[ORGID] is the identifier for the organism database, e.g. ECOLI, META, AFER243159
[2|3] specifies whether data should use BioPAX Level 2 or Level 3. If the type argument is omitted, BioPAX Level 3 will be generated.
[PATHWAY] is in the BioCyc identifier for the pathway, e.g. GLYCOLYSIS, ARGSYN-PWY, PWY0-1299
Example URLs:

https://websvc.biocyc.org/BSUB/pathway-biopax?type=3&object=GLYCOLYSIS
Retrieve the glycolysis pathway in Bacillus subtilis in BioPAX Level 3 format.
https://websvc.biocyc.org/AFER243159/pathway-biopax?type=2&object=CYSTSYN-PWY
Retrieve the cysteine biosynthesis pathway in Acidithiobacillus ferrooxidans in BioPAX Level 2 format.
https://websvc.biocyc.org/META/pathway-biopax?object=PWY-5025
Retrieve the IAA biosynthesis pathway in MetaCyc in BioPAX Level 3 format.
Search and Retrieve a Set of Objects
Retrieving a Set of Objects Using the Pathway Tools API Functions
The set of Pathway Tools API functions were defined to allow users who have downloaded and installed the Pathway Tools software locally (why do this?) to write programs that operate on the data. A subset of these API functions have been made available via the web services interface.

The URL to issue an API query that returns a list of objects in ptools-xml format is

https://websvc.biocyc.org/apixml?fn=[API-FUNCTION]&id=[ORGID]:[OBJECT-ID]&detail=[none|low|full]

where [ORGID], [OBJECT-ID] and detail level are as for the single object queries described above, and [API-FUNCTION] is one of the following:

all-products-of-gene
binding-site-transcription-factors
chromosome-of-gene
compounds-of-pathway
containers-of
containing-tus
direct-activators
direct-inhibitors
enzymes-of-gene
enzymes-of-pathway
enzymes-of-reaction
genes-of-pathway
genes-of-protein
genes-of-reaction
genes-regulated-by-gene
genes-regulating-gene
get-class-all-instances
get-class-all-subs
get-class-direct-supers
get-class-direct-subs
modified-containers
modified-forms
monomers-of-protein
pathways-of-compound
pathways-of-gene
reactions-of-compound
reactions-of-enzyme
reactions-of-gene
regulator-proteins-of-transcription-unit
regulon-of-protein
substrates-of-reaction
top-containers
transcription-unit-activators
transcription-unit-binding-sites
transcription-unit-genes
transcription-unit-inhibitors
transcription-unit-mrna-binding-sites
transcription-unit-promoter
transcription-unit-terminators
transcription-unit-transcription-factors
transcription-units-of-gene
transcription-units-of-protein
A more detailed description of each API function is available here.

Example URLs:

https://websvc.biocyc.org/apixml?fn=genes-of-pathway&id=BSUB:GLYCOLYSIS
Retrieve the set of genes that participate in the glycolysis pathway in Bacillus subtilis.
https://websvc.biocyc.org/apixml?fn=genes-regulated-by-gene&id=ECOLI:EG10164&detail=none
Retrieve the set of genes (IDs only) regulated by the crp gene in EcoCyc.
https://websvc.biocyc.org/apixml?fn=enzymes-of-reaction&id=META:TRYPSYN-RXN&detail=full
Get detailed information on all enzymes in MetaCyc that catalyze the tryptophan synthase reaction.
https://websvc.biocyc.org/apixml?fn=get-class-all-instances&id=META:Genes&detail=full
Get detailed information on all genes in MetaCyc.
Retrieving the Set of Objects Returned by a BioVelo Query
The URL to issue a BioVelo query that returns a list of objects in ptools-xml format is

https://websvc.biocyc.org/xmlquery?[QUERY] or
https://websvc.biocyc.org/xmlquery?query=[QUERY]&detail=[none|low|full]

where

[QUERY] is a properly escaped BioVelo query string that returns a single list of Pathway Tools objects (it is possible to create BioVelo queries that return multi-column tables -- these are not appropriate as input for this web service). For more information about constructing BioVelo queries, see the Guide to the BioVelo Query Language.
[none|low|full] indicates whether the returned output should contain no detail (i.e. only object identifiers and links will be included), low detail (names and a handful of other attributes will be included) or full detail for the matching objects. If no detail parameter is supplied, the request defaults to low detail.
NOTE: the BioVelo query needs to be URL encoded if it contains any non-alphanumeric characters. You might look here, or here if you are unfamiliar with URL encoding.

Example URLs:

https://websvc.biocyc.org/xmlquery?%5bx%3ax%3c-ecoli%5e%5epathways%5d
Retrieve the complete set of pathways in the EcoCyc database at low detail, where the query [x:x<-ecoli^^pathways]
has been URL encoded as %5bx%3ax%3c-ecoli%5e%5epathways%5d.
https://websvc.biocyc.org/xmlquery?query=[x:x<-ecoli^^genes,x^name%3D"trpA"]&detail=full
Retrieve the gene (or genes) in the EcoCyc database with the name "trpA" in full detail,
where the query [x:x<-ecoli^^genes,x^name%3D"trpA"]
has been encoded as %5bx%3ax%3c-ecoli%5e%5egenes%2cx%5ename%253D%22trpA%22%5d.
https://websvc.biocyc.org/xmlquery?dbs
Retrieve the set of available organism databases, and the query requires no URL encoding.
https://websvc.biocyc.org/xmlquery?%5bx%3ay%3a%3dbsub~argsynbsub-pwy%2cx%3c-(enzymes-of-pathway+y)%5d
Retrieve the set of enzymes that participate in the arginine biosynthesis pathway in Bacillus subtilis,
where the query [x:y:=bsub~argsynbsub-pwy,x<-(enzymes-of-pathway y)]
has been URL encoded as %5bx%3ay%3a%3dbsub~argsynbsub-pwy%2cx%3c-(enzymes-of-pathway+y)%5d.
https://websvc.biocyc.org/xmlquery?%5bx%3ax%3c-meta%5e%5eproteins%2c+%22aspartate%22+instringci+x%5enames%26%22kinase%22+instringci+x%5enames%5d
Retrieve the set of proteins in MetaCyc that have the words "aspartate" and "kinase" in their common-name or synonyms
where the query [x:x<-meta^^proteins, "aspartate" instringci x^names&"kinase" instringci x^names]
has been URL encoded as %5bx%3ax%3c-meta%5e%5eproteins%2c+%22aspartate%22+instringci+x%5enames%26%22kinase%22+instringci+x%5enames%5d.
https://websvc.biocyc.org/xmlquery?query=%5bx%3ax%3c-ecoli%5e%5epathways%2cecoli~EG10258+in+(pathway-to-genes+x)%5d&detail=none
Retrieve the identifiers of all pathways containing the eno (enolase) gene in Escherichia coli,
where the query [x:x<-ecoli^^pathways,ecoli~EG10258 in (pathway-to-genes x)]
has been encoded as %5bx%3ax%3c-ecoli%5e%5epathways%2cecoli~EG10258+in+(pathway-to-genes+x)%5d.
Retrieve BioCyc Objects from a Class by Name or Accession
Objects from a specified class in a PGDB can be retrieved by name or accession:
https://websvc.biocyc.org/[ORGID]/name-search?object=[NAME]&class=[CLASS] or
https://websvc.biocyc.org/[ORGID]/name-search?object=[NAME]&class=[CLASS]&fmt=[json|xml]

where

[ORGID] is the identifier for the organism database, e.g. ECOLI, META, AFER243159.
[NAME] is the name or accession of the object to retrieve.
[CLASS] is the class from which the object is to be retrieved, e.g. Genes.
fmt=json requests output in the JSON format; default output is in tab-delimited format.
fmt=xml requests output in ptools-XML format; the output will include the full description of the matching object.
This web service can be used if the user is retrieving the results by either the "POST" or "GET" method.
Example URLs:

https://websvc.biocyc.org//ECOLI/name-search?object=trpa&class=Genes
Retrieve Genes called trpA in EcoCyc.
https://websvc.biocyc.org/ECOLI/name-search?object=trpa&class=Genes&fmt=json
Retrieve Genes called trpA in EcoCyc in JSON format.
https://websvc.biocyc.org/ECOLI/name-search?object=b1234&class=Genes&fmt=xml
Retrieve Genes with accession b1234 in EcoCyc in ptools-XML format.
The default output is a tab-delimited format as follows. If no results are found then the string EMPTY-RESULT is returned.

EG11024	trpA
where

column 1 : BioCyc Object ID
column 2 : Common Name

JSON format:

{"RESULTS":[{"OBJECT-ID":"EG11024","COMMON-NAME":"trpA"}]}

Retrieving a BioCyc Object Given a Foreign ID
This service finds the BioCyc ID of an object given a foreign ID, that is, an identifier of that object in an external database. This service depends on the foreign ID being stored in the DB-Links slot of a BioCyc object. All DB-Links in a given PGDB are searched. Example foreign IDs: Uniprot IDs (P11509) for proteins; KEGG IDs (C00033) for compounds, etc.; Rhea IDs (10743) for reactions.
BioCyc objects such as pathways, reactions, compounds, genes, proteins, RNAs, or transcription-units in a PGDB can be retrieved given a foreign identifier that matches that object.

The URLs to search for objects based on foreign identifier are

https://websvc.biocyc.org/[ORGID]/foreignid?ids=[DATABASE-NAME]:[FOREIGNID] or
https://websvc.biocyc.org/[ORGID]/foreignid?ids=[DATABASE-NAME]:[FOREIGNID]&fmt=[json|xml]

where

[ORGID] is the identifier for the organism database, e.g. ECOLI, META, AFER243159.
[DATABASE-NAME] is the name of external database, e.g. KEGG, PUBCHEM, UNIPROT.
[FOREIGNID] is the external database identifier. The foreignid input may contain one or more values, which may be separated by commas. Note that foreignids are case-sensitive.
fmt=json requests output in the JSON format; default output is in tab-delimited format.
fmt=xml requests output in ptools-XML format; the output will include the full description of all matching objects (including all links to external databases), but will not indicate which of the input ids match which objects, and will not identify any input ids that did not match any objects.
This web service can be used if the user is retrieving the results by either the "POST" or "GET" method.
Example URLs:

https://websvc.biocyc.org/META/foreignid?ids=Kegg:C00849
Retrieve the object that contains a KEGG link of C00849 in MetaCyc.
https://websvc.biocyc.org/ECOLI/foreignid?ids=Kegg:C00849,EcoGene:EG11025
Retrieve objects in EcoCyc that contain a KEGG link of C00849 or an EcoGene link of EG11025.
https://websvc.biocyc.org/BSUB/foreignid?ids=UniProt:P05653&fmt=json
Retrieve the object that contains a UniProt link of P05653 in Bacillus subtilis in JSON format.
https://websvc.biocyc.org/BSUB/foreignid?ids=UniProt:P05653&fmt=xml
Retrieve the object that contains a UniProt link of P05653 in Bacillus subtilis in ptools-XML format.
The default output is a tab-delimited format as follows:

Kegg:C00849	1	ETHYLACETATE
EcoGene:EG11025	1	EG11025
where

column 1 : Input Foreign ID.
column 2 : 1 means that a valid object that matches the foreign ID was found and 0 means that no object that matches the foreign ID was found.
column 3 : The BioCyc identifier of the object that matches the foreign ID.

JSON format:

[{"INPUT":"UniProt:P05653","STATUS":1,"RESULTS":[{"ID":"BSU00070-MONOMER"}]}]


Retrieving a BioCyc Object Given a Chemical Formula
This service finds the BioCyc IDs of all metabolites that exactly match a supplied chemical formula. This service consults the Chemical-Formula slot of BioCyc compounds.
The URLs to search metabolites based on chemical formula are

https://websvc.biocyc.org/[ORGID]/CF?cfs=[CHEMICAL-FORMULA] or
https://websvc.biocyc.org/[ORGID]/CF?cfs=[CHEMICAL-FORMULA]&fmt=[json|xml]

where

[ORGID] is the identifier for the organism database, e.g. ECOLI, META, AFER243159.
[CHEMICAL-FORMULA] is the chemical formula to be matched, in the format ([element-symbol][coefficient])+. The chemical formula input may contain one or more values, which may be separated by commas. The search will return all metabolites in the specified database that have the exact chemical formula(s) provided here.
fmt=json requests output in JSON format; default output is in tab-delimited format
fmt=xml requests output in ptools-XML format; the output will include the full description of all matching compounds (including chemical formulas), but will not indicate which of the input formulas (if multiple are supplied) match which compounds, and will not identify any input formulas that did not match any compounds.
This web service can be used if the user is retrieving the results by either the "POST" or "GET" method. If multiple Chemical formulas are provided than this web service will search for compounds matching each of the chemical formula.
Example URLs:

https://websvc.biocyc.org/META/CF?cfs=C6H6
Retrieve the compound that contains a chemical formula of C6H6 in MetaCyc.
https://websvc.biocyc.org/ECOLI/CF?cfs=C6H12O6,H2O
Retrieve the compounds that contain the chemical formulas C6H12O6 and H2O in EcoCyc.
https://websvc.biocyc.org/BSUB/CF?cfs=C6H12O6&fmt=json
Retrieve the compounds that contains a Chemical formula of C6H12O6 in Bacillus subtilis in JSON format.
https://websvc.biocyc.org/BSUB/CF?cfs=C6H12O6&fmt=xml
Retrieve the compounds that contains a Chemical formula of C6H12O6 in Bacillus subtilis in ptools-XML format.
The default output is a tab-delimited format as follows:

C6H6	1	BENZENE		benzene
H2O3	0
where

column 1 : Input Chemical Formula.
column 2 : 1 means that a valid compound that matches the chemical formula was found; 0 means that no match was found.
column 3 : The BioCyc identifier of the compound that matches the chemical formula.
column 4 : The common-name of the compound.

JSON format:

[{"INPUT":"H2O2","STATUS":1,"RESULTS":[{"ID":"HYDROGEN-PEROXIDE","NAME":"hydrogen peroxide"}]},{"INPUT":"H2O3","STATUS":0}]


Retrieving a BioCyc Compound given a Monoisotopic Weight and Tolerance
Compounds in a PGDB can be retrieved given a monoisotopic molecular weight and tolerance.
The URLs to search for compounds based on monoisotopic molecular weight and tolerance are

https://websvc.biocyc.org/[ORGID]/monoisotopicwt?wts=[MONOISOTOPICMW]&tol=[TOLERANCE] or
https://websvc.biocyc.org/[ORGID]/monoisotopicwt?wts=[MONOISOTOPICMW]&tol=[TOLERANCE]&fmt=[json|xml]

where

[ORGID] is the identifier for the organism database, e.g. ECOLI, META, AFER243159.
[MONOISOTOPICMW] is a floating point number, standing for the monoisotopic molecular weight in daltons. One or more monoisotopic molecular weights can be supplied, where multiple values are separated by commas.
[TOLERANCE] is the search tolerance +/- in ppm.
fmt=json requests output in the JSON format; default output is in tab-delimited format
fmt=xml requests output in ptools-XML format; the output will include the full description of all matching compounds (including their monoisotopic molecular weights), but will not indicate which of the input weights (if multiple are supplied) match which compounds, and will not identify any input weights that did not match any compounds.
This web service can be used if the user is retrieving the results by either the "POST" or "GET" method.
Example URLs:

https://websvc.biocyc.org//META/monoisotopicwt?wts=123.009&tol=5
Retrieve the compounds that have a monoisotopic molecular weight of 123.009 and tolerance of 5 in MetaCyc.
https://websvc.biocyc.org/ECOLI/monoisotopicwt?wts=240.063,88.052&tol=5
Retrieve the compounds that have a monoisotopic molecular weight of either 240.063 or 88.052 and tolerance of 5 in EcoCyc.
https://websvc.biocyc.org/BSUB/monoisotopicwt?wts=169.988&tol=5&fmt=JSON
Retrieve the compounds that have a monoisotopic molecular weight of 169.988 and tolerance of 5 in Bacillus subtilis in JSON format.
https://websvc.biocyc.org/BSUB/monoisotopicwt?wts=169.988&tol=5&fmt=xml
Retrieve the compounds that have a monoisotopic molecular weight of 169.988 and tolerance of 5 in Bacillus subtilis in ptools-XML format.
The default output is a tab-delimited format as follows:

123.009	1	123.008705	3-chloro-L-alanine	CHLORALAN-CPD
123.009	1	123.008705	3-chloro-D-alanine	3-CHLORO-D-ALANINE
123.009	1	123.008705	2-chloro-L-alanine	CPD0-1475
123.009	1	123.008705	3-chloro-DL-alanine	3-CHLORO-DL-ALANINE
56	0			
where

column 1 : Input monoisotopic molecular weight for this search.
column 2 : 1 means the query was successful and 0 means the query found no match.
column 3 : Monoisotopic molecular weight of the compound that is stored in the PGDB
column 4 : Compound name
column 5 : The BioCyc identifier of the compound.

JSON format:

[{"INPUT":"123.009","STATUS":1,"RESULTS":[{"MW":123.008705,"NAME":"3-chloro-L-alanine","ID":"CHLORALAN-CPD"},{"MW":123.008705,"NAME":"3-chloro-D-alanine","ID":"3-CHLORO-D-ALANINE"},{"MW":123.008705,"NAME":"2-chloro-L-alanine","ID":"CPD0-1475"},{"MW":123.008705,"NAME":"3-chloro-DL-alanine","ID":"3-CHLORO-DL-ALANINE"}]},{"INPUT":"56","STATUS":0}]

Metabolite Translation Service
This web service translates metabolite names, identifiers, InChI strings, InChI keys, monoisotopic molecular weights, and molecular formula, between metabolite databases. Its input is a set of lines, one line per metabolite. Each line of the file contains one or more metabolite names, identifiers, and an optional InChI string, InChI key, monoisotopic molecular weight, and chemical formula. The service looks up each of the preceding fields within the specified BioCyc database. Three cases are possible for each line:

1) None of the provided fields is recognized in the specified BioCyc database, in which case "unknown" is returned, along with the unknown input fields.
2) All of the recognized names, identifiers, InChI string, InChI key, monoisotopic weight, and chemical formula match a single metabolite, in which case "successful" is returned, along with the following tab-separated fields:

BioCyc ID of matching metabolite
BioCyc common name of matching metabolite
Additional identifiers present in BioCyc from other databases for the matching metabolite
3) The recognized names, identifiers, InChI string, InChI key, monoisotopic weight, and chemical formula match more than one metabolite, in which case "ambiguous" is returned, along with the ambiguous fields from the input.

The default input for the file or for pasting data is a tab-delimited format.

Example:

Kegg:C00001	PubChem:125
TRP
ACET
The default output is a tab-delimited format as follows:

ambiguous	BioCyc:WATER		BioCyc:4-HYDROXY-BENZYL-ALCOHOL	
success		BioCyc:TRP	L-tryptophan	MetaboLights:MTBLC57912	HMDB:HMDB00929	IAF1260:33772	ChEBI:57912	PubChem:6923516	KEGG:C00078	CAS:73-22-3
success		BioCyc:ACET	acetate		MetaboLights:MTBLC30089	HMDB:HMDB00042	DrugBank:DB03166	IAF1260:33590	ChemSpider:170	PubChem:175	ChEBI:30089	CAS:64-19-7	KEGG:C00033	CAS:71-50-1


To access the metabolite translation service interactively, click here


Retrieve BioCyc Objects by SMILES String
Retrieve compounds whose chemical structure matches a given SMILES string:
https://websvc.biocyc.org/[ORGID]/smiles-search?smiles=[SMILES]&exact=[EXACT] or
https://websvc.biocyc.org/[ORGID]/smiles-search?smiles=[SMILES]&exact=[EXACT]&fmt=json

where

[ORGID] is the identifier for the organism database, e.g. ECOLI, META, AFER243159.
[SMILES] is the SMILES string, e.g. CN, CO.
[EXACT] is a Boolean, either T (true) or F (false).
fmt=json requests output in the JSON format; default output is in tab-delimited format.
This web service can be used if the user is retrieving the results by either the "POST" or "GET" method.
Example URLs:

https://websvc.biocyc.org/ECO/smiles-search?smiles=CO&exact=F
Retrieve compounds containing CO, inexact match.
https://websvc.biocyc.org/ECO/smiles-search?smiles=CO&exact=F&fmt=json
Retrieve compounds containing CO, inexact match, using JSON format.
The default output is a tab-delimited format as follows. If no results are found then the string EMPTY-RESULT is returned.

CPD-22525	p-aminophenyl-β-D-glucoside	        C(O)[C@@H]2([C@@H](O)[C@H](O)[C@@H](O)[C@H](OC1(\C=C/C(/N)=C\C=1))O2)
CPD-8122	molybdopterin adenine dinucleotide	C(OP(=O)([O-])OP(OC[C@H]2(O[C@H]1(NC3(/N=C(N)\NC(=O)C(/N[C@H]1C(\S)=C2\[S-])=3))))(=O)[O-])[C@H]6(O[C@@H](N4(C5(\C(\N=C/4)=C(N)/N=C\N=5)))[C@H](O)[C@H](O)6)
... 
where

column 1 : BioCyc Object ID
column 2 : Common Name
column 3 : Complete SMILES string

JSON format:

{"RESULTS":[{"OBJECT-ID":"CPD-22525","COMMON-NAME":"p-aminophenyl-β-D-glucoside","SMILES":"C(O)[C@@H]2([C@@H](O)[C@H](O)[C@@H](O)[C@H](OC1(\\C=C\/C(\/N)=C\\C=1))O2)"},{"OBJECT-ID":"CPD-8122","COMMON-NAME":"molybdopterin adenine dinucleotide","SMILES":"C(OP(=O)([O-])OP(OC[C@H]2(O[C@H]1(NC3(\/N=C(N)\\NC(=O)C(\/N[C@H]1C(\\S)=C2\\[S-])=3))))(=O)[O-])[C@H]6(O[C@@H](N4(C5(\\C(\\N=C\/4)=C(N)\/N=C\\N=5)))[C@H](O)[C@H](O)6)"},...]}

Retrieve BioCyc Objects by InChI String
Retrieve chemical compounds that match a given InChI string:
https://websvc.biocyc.org/[ORGID]/inchi-search?inchi=[INCHI]&exact=[EXACT] or
https://websvc.biocyc.org/[ORGID]/inchi-search?inchi=[INCHI]&exact=[EXACT]&fmt=json

where

[ORGID] is the identifier for the organism database, e.g. ECOLI, META, AFER243159.
[INCHI] is the InChI string, e.g. 1S/C4H5N3O/c5-3-1-2-6-4(8)7-3/h1-2H,(H3,5,6,7,8)
fmt=json requests output in the JSON format; default output is in tab-delimited format.
This web service can be used if the user is retrieving the results by either the "POST" or "GET" method.
Example URLs:

https://websvc.biocyc.org/ECO/inchi-search?inchi=1S/C4H5N3O/c5-3-1-2-6-4(8)7-3/h1-2H,(H3,5,6,7,8)
Retrieve compounds with INCHIstring 1S/C4H5N3O/c5-3-1-2-6-4(8)7-3/h1-2H,(H3,5,6,7,8).
https://websvc.biocyc.org/ECO/inchi-search?inchi=1S/C4H5N3O/c5-3-1-2-6-4(8)7-3/h1-2H,(H3,5,6,7,8)&fmt=json
Retrieve compounds with INCHIstring 1S/C4H5N3O/c5-3-1-2-6-4(8)7-3/h1-2H,(H3,5,6,7,8), using JSON format.
The default output is a tab-delimited format as follows. If no results are found then the string EMPTY-RESULT is returned.

CYTOSINE	cytosine	InChI=1S/C4H5N3O/c5-3-1-2-6-4(8)7-3/h1-2H,(H3,5,6,7,8)	(1S/C4H5N3O/c5-3-1-2-6-4(8)7-3/h1-2H,(H3,5,6,7,8) 1S/C4H5N3O/c5-3-1-2-6-4(8)7-3/h1 1S/C4H5N3O/c5)
where

column 1 : BioCyc Object ID
column 2 : Common Name
column 3 : INCHI String
column 4 : List of INCHI Keys

JSON format:

{"RESULTS":[{"OBJECT-ID":"CYTOSINE","COMMON-NAME":"cytosine","INCHI-STRING":"InChI=1S\/C4H5N3O\/c5-3-1-2-6-4(8)7-3\/h1-2H,(H3,5,6,7,8)","INCHI-KEYS":["1S\/C4H5N3O\/c5-3-1-2-6-4(8)7-3\/h1-2H,(H3,5,6,7,8)","1S\/C4H5N3O\/c5-3-1-2-6-4(8)7-3\/h1","1S\/C4H5N3O\/c5"]}]}

Visualization Services
Generate Compound Image
Compound images may be retrieved in GIF format for incorporation into other pages. The basic URL to retrieve a compound image is:
https://websvc.biocyc.org/[ORGID]/diagram-only?type=COMPOUND&object=[COMPOUND]

where

[ORGID] is the identifier for the organism database, e.g. ECOLI, META, AFER243159
[COMPOUND] is in the BioCyc identifier for the compound, e.g. TRP, CPD-560
Example URL:

https://websvc.biocyc.org/ECOLI/diagram-only?type=COMPOUND&object=TRP
Generate Pathway Diagram
Pathway images may be generated in GIF format for incorporation into other pages. These images are fully customizable as to the level of detail and which elements are included in the diagrams. URLs to generate a pathway image are of the form:
https://websvc.biocyc.org/[ORGID]/diagram-only?type=PATHWAY&object=[PATHWAY]&[PARAMETER]=[VALUE]

where

[ORGID] is the identifier for the organism database, e.g. ECOLI, META, AFER243159
[PATHWAY] is in the BioCyc identifier for the pathway, e.g. GLYCOLYSIS, ARGSYN-PWY, PWY0-1299
Multiple parameter/value pairs can be specified
The following parameters provide customization options for the pathway diagrams:


Parameters to Customize Pathway Appearance
Parameter Name	Possible Values	Description	Default
detail-level	0,1,2,3,4	Specifies the general detail level for the pathway. The specific elements shown or hidden can be individually overridden using the other parameters below.
0 - Minimal detail, only start/end/branchpoint metabolites shown
1 - All main compounds (those along the main pathway backbone)
2 - All main and side compounds, enzyme and gene names, and EC#s
3 - Compound structures for most main compounds
4 - Compound structures for most main and side compounds	0, 1 or 2, depending on the size and complexity of the individual pathway.
enz	y, n	Should enzyme names be shown?	Depends on detail level
gene	y, n	Should gene names be shown?	Depends on detail level
ec	y, n	Should EC numbers be shown?	Depends on detail level
secs	y, n	Should side compounds be shown?	Depends on detail level
mstruct	none, most, all	Should compound structures for main compounds be shown? When most is selected, structures for some very common compounds, such as ATP, NAD, etc. are omitted.	none or most, depending on detail level
sstruct	none, most, all	Should compound structures for side compounds be shown? When most is selected, structures for some very common compounds, such as ATP, NAD, etc. are omitted.	none or most, depending on detail level
reglinks	y, n	If n, omit enzyme regulation icons and feedback inhibition links	y
nolinks	y, n	If y, suppress showing links to other pathways	n
pfontsize	tiny, very-small, small, normal, large, very-large	Font size used in diagram	very-small
bgcolor	w, cb, g, bw, tr	The color scheme for the diagram:
w: Colors on a white background
cb: Colors on a black background
g: Colors on a gray background
bw: Black on a white background
tr: Same colors as for w, but on a transparent background	tr
linear	snake, horizontal, vertical	Specify whether linear pathways are drawn horizontally, vertically, or in a horizontal back-and-forth "snake"-like fashion. Note that only pathways that contain no cycles or branchpoints (including branchpoints resulting from links to/from other pathways) are recognized as linear pathways and therefore sensitive to this parameter.	snake
Example URLs:

https://websvc.biocyc.org/ECOLI/diagram-only?type=PATHWAY&object=GLYCOLYSIS
https://websvc.biocyc.org/ECOLI/diagram-only?type=PATHWAY&object=GLYCOLYSIS&detail-level=1&reglinks=n&nolinks=y&bgcolor=bw
https://websvc.biocyc.org/BSUB/diagram-only?type=PATHWAY&object=TYRSYN&enz=n&gene=y&ec=n&mstruct=most


Genome Browser Visualization
Genome Explorer can be used to explore one replicon (chromosome or plasmid) at a time, to align conserved regions of multiple replicons, and to analyze high-throughput datasets using tracks. The following are a set of tools that provide customization options for the genome browser via url parameters to achieve the desired effect.


Single Genome Range Request:

Display a specified base-pair range for a specified replicon.

Parameters:
https://websvc.biocyc.org/genbro/genbro.shtml?orgid=[ORGID]&replicon=[REPLICON]&bp-range=[BP-RANGE]

OrgID: The organism's ID.
Example: ECOLI

Replicon: The replicon's ID.
Example: COLI-K12

Bp-Range: A coordinate range (e.g., START/END) to display that region.
Example: 1000/20000

Example URL:
https://websvc.biocyc.org/genbro/genbro.shtml?orgid=ECOLI&replicon=COLI-K12&bp-range=500000/1000000


Specific Gene in Genome Request:
Display a region of a specified replicon centered on a specified gene.

Parameters:
https://websvc.biocyc.org/genbro/genbro.shtml?orgid=[ORGID]&replicon=[REPLICON]&gene=[GENE]

OrgID: The organism's ID.
Example: ECOLI

Replicon: The replicon's ID.
Example: COLI-K12

Gene: The gene's ID.
Example: EG10240

Example URL: https://websvc.biocyc.org/genbro/genbro.shtml?orgid=ECOLI&replicon=COLI-K12&gene=EG10240


Comparative Genome Explorer Request:
Produce a comparative view of several specified replicons aligned at orthologs of a specified gene.

Parameters: https://websvc.biocyc.org/genbro/genbro.shtml?lead-orgid=[LEAD-ORGID]&lead-genes=[LEADGENE]&replicon=[REPLICON]&orgids=[ORGIDS]

Lead OrgID: The leading organism's ID
Example: ECOLI

Lead Gene: The lead gene's ID
Example: EG10240

Replicon: The replicon's ID
Example: COLI-K12

OrgIDs: A list or organism ids, should be seperated by commas. (e.g., ORGID,ORGID2,ORGID3)
Example: ECOLI,GCF_000340255

Example URL: https://websvc.biocyc.org/genbro/ortho.shtml?lead-orgid=ECOLI&lead-genes=EG10240&replicon=COLI-K12&orgids=ECOLI,GCF_000340255





Omics Visualization Services
Pathway Diagram Painted with Omics Data
Pathway pages can be displayed with omics data superimposed on the pathway diagrams. This operation navigates to the web page for the requested pathway. URLs to generate a pathway web page with overlayed omics data are of the form:
https://websvc.biocyc.org/[ORGID]/new-image?type=PATHWAY&object=[PATHWAY]&[PARAMETER]=[VALUE]

where

[ORGID] is the identifier for the organism database, e.g. ECOLI, BSUB, AFER243159
[PATHWAY] is in the BioCyc identifier for the pathway, e.g. GLYCOLYSIS, ARGSYN-PWY, PWY0-1299
Multiple parameter/value pairs can be specified, but two are required in order to show omics data:
url (for GET) or datafile (for POST) -- The URL of the omics data file
column1 -- An integer specifying the column number for the data from the omics data file that is to be displayed (column zero is the first column)
The table below outlines the possible paramaters that control the mapping of omics data to pathway diagrams, as well as on the cellular overview diagram (see the next section). In addition, the options detailed in the previous table can be used to customize the rendering of the pathway itself.

Users can submit omics data using a POST request (which allows the omics data to be private). For our examples below we make a data file available on a Pathway Tools web server.


Parameters to Customize Pathway and Cellular Overview Omics-Data Visualizations
Parameter Name	Possible Values	Pathway (P), Table (T) or Cellular Overview (O)	Description	Default
datafile	A filename	P, T, O	This parameter is REQUIRED for a POST, and should not be used for a GET. Its value should be a filename.
The file should be tab-delimited, the first column (often referred to as Column 0) being the id of the gene/protein/compound/reaction/other to which the data in the remaining columns applies.	REQUIRED
for POST
no default
url	An example file supplied by the Pathway Tools website administrator	P, T, O	This parameter is REQUIRED for a GET, and should not be used for a POST. Its value should name an example data file supplied by the Pathway Tools website administrator.	REQUIRED
for GET
no default
class	gene, protein, compound, reaction, 'NIL'	P, T, O	gene Gene names and/or identifiers
protein Protein names and/or identifiers
compound Compound names and/or identifiers
reaction Reaction identifiers and/or EC numbers
'NIL' Any of the above

The given value declares the type(s) of the ids found in column 0 (the first column) of the data file.	gene
column1	numbers separated by space or comma AND/OR,
range of numbers (number pair separated by hyphen)	P, T, O	data columns from datafile to be painted on pathway graphic (e.g. &column1=4-7 9 would select data from columns 4, 5, 6, 7, and 9 in each row of the file.)	REQUIRED
no default
expressiontype	relative, absolute	P, T, O	relative: the numerical data represents a ratio (relative to some other dataset or to a control) and is therefore centered around either 0 or 1 (see log parameter below).
absolute: the numerical data represents absolute values, e.g. intensities or concentrations, and are all non-negative.	relative
log	on, off	P, T, O	on: 0-centered data
off: 1-centered data
0-centered scale: implies that the numerical data of your file can contain positive and negative values. The value 0 is considered to be the center of the numerical values provided in your data file.
1-centered scale: implies that any negative or zero values in your data file should be skipped. Moreover, the data is centered around the value 1 using a log scale. For example, the value 0.1 is considered to be at the same distance to 1 as the value 10. So, a logarithm of base 10 is applied to your data before the linear coloring mapping is applied.
on
color	default, specify, rbg, rbg-cutoff, 3-color	P, T	
default : Orange to gray to blue from data
specify : Orange to gray to blue with a maximum cutoff (maxcutoff parameter detailed in this table)
rbg : Red to green from data
rbg-cutoff : Red to green with a maximum cutoff (maxcutoff parameter detailed in this table)
3-color : Three-color display with a threshold
Data values are divided into color bins. You can choose between a color scheme that ranges from orange (most positive values) through gray in the center to blue (most negative values) or from red through blue in the center to green and yellow. For each of these color schemes, two options are available:
The color bins range over the entire spectrum, and the cutoff values for the color bins are derived from the data itself. This means that different experiments could be displayed using different color schemes, making it difficult to directly compare them.
You may specify a value for the maximum value cutoff (maxcutoff parameter) bin. All displays that use the same maximum value cutoff will use the same color scheme (assuming other settings are the same), and are therefore directly comparable. The maximum cutoff value should be a number, e.g. 2 or 10, etc. All data values greater than the maximum cutoff value will be displayed in the highest bin color.
A final alternative is to use only three color bins, red for data values that exceed some threshold (see parameter below), purple for data values less than the inverse of that threshold, and gray for values in between. The threshold value should be a number, e.g. 2 or 10.	default
threshold	a numeric value	P, T	Required if 3-color specified. See discussion in 'color' parameter section of this table.	no default
maxcutoff	a numeric value	P, T	Required if specify or rbg-cutoff specified. See discussion in 'color' parameter section of this table.	no default
omicsPopups	on, unspecified	P	If present, display the omics data in popup windows, usually one popup per reaction or metabolite in the pathway. This is the only way to show data in multiple columns. A single column of data can either be shown in popups or (if unspecified) as color coded buttons on the diagram.	unspecified
defaultPopup	bar, plot, heat	P	The style of data display in the omics popups: a bar graph, an x-y plot, or a heat map.	bar

A single column of gene data from the file displayed with colored "buttons":
https://websvc.biocyc.org/ECOLI/new-image?type=PATHWAY&object=TCA&detail-level=2&url=/expr-examples/sample.dat&expressiontype=relative&numcolumns=1&column1=18&log=on&class=gene&color=Default

A single column of gene data from the file displayed with popups:
https://websvc.biocyc.org/ECOLI/new-image?type=PATHWAY&object=TCA&detail-level=2&url=expr-examples/sample.dat&expressiontype=relative&numcolumns=1&column1=18&log=on&class=gene&color=Default&omicsPopups=on

Five columns of compound (metabolomics) data from the file displayed with popups and heat map as default popup style:
https://websvc.biocyc.org/HUMAN/new-image?type=PATHWAY&object=PWY66-422&url=expr-examples/human-metabolomics-syn.txt&expressiontype=relative&column1=1-5&log=on&class=compound&color=Default&omicsPopups=on&defaultPopup=heat
Table of Pathway Diagrams Painted with Omics Data
This operation navigates to a web page containing a table of low-detail pathway diagrams, colored to show omics data values. To submit data via the HTTP GET method, use the following URL format:
https://websvc.biocyc.org/[ORGID]/overview-expression-map?pathways=[PATHWAYS]&[PARAMETER]=[VALUE]

where

[ORGID] is the identifier for the organism database, e.g. ECOLI, BSUB, AFER243159
[PATHWAYS] is a comma-separated list of BioCyc identifiers for the desired pathways, e.g. GLYCOLYSIS, ARGSYN-PWY, PWY0-1299
Multiple parameter/value pairs can be specified. If none are specified, then the table will be generated without omics data.
The list of available parameters are those marked T in the above Omics Data Parameters Table. The pathway diagrams cannot be customized, and any pathway customization parameters will be ignored. This view is more useful for gene expression data than for metabolomics data, as any side metabolites will not be visible on the pathway diagrams. Popups are also not available with this view -- if multiple data columns are specified, the resulting table will contain a column of pathway diagrams for each data column. The same URL and parameters are used for GET and POST queries except that GET queries must supply the url parameter and POST queries must supply the datafile parameter to upload data from a local file.

Example URL:

https://websvc.biocyc.org/ECOLI/overview-expression-map?url=expr-examples/sample.dat&expressiontype=relative&column1=18&log=on&class=gene&pathways=GLYCOLYSIS,TCA,TRPSYN-PWY
Cellular Overview Visualization Services
The following two services navigate to a web page containing the cellular overview (metabolic map) diagram with specified data overlays. For the highlighting service the objects to be highlighted are specified in the URL itself, and highlighting colors are optionally provided. For the omics painting service the objects to be highlighted are specified in a POST, and data overlays are derived from omics data values that are mapped to a color scale.
Highlight/Paint Objects on Metabolic Network Diagram
Objects can be specified for highlighting by name, ID or substring. URLs to highlight specified objects are of the form:
https://websvc.biocyc.org/overviewsWeb/celOv.shtml?orgid=[ORGID]&zoomlevel=[0-3]&[OP]=[VALUES]

where

[ORGID] is the identifier for the organism database, e.g. ECOLI, BSUB, AFER243159
zoomlevel is an integer that indicates the extent to which the diagram has been zoomed. Zoomlevel 0 is the lowest overview magnification, level 3 is the highest magnification. If not supplied, the zoomlevel defaults to 0.
[OP] specifies a highlighting operation
[VALUES] consists of one or more values, where multiple values are separated by URL-encoded spaces (i.e., each space is represented by the characters %20). Each value is of the form [LABEL] [NUMBER] where [LABEL] is the name, ID or substring indicating the object to be highlighted, and [NUMBER] is an optional omics data value that is mapped to a highlighting color ([LABEL] and [NUMBER] are separated by a URL-encoded space). An example [OP] plus [VALUES] is "xnids=b1234 35 b3320 86 b2310 45" where b1234 is a gene identifier and 35 is an associated data value.
Multiple highlighting operations can be specified in a single URL. The possible highlighting operations are described in the following table. Note that all these operations correspond to the operations available from the top menu bar when a Cellular Overview diagram is displayed. The operation xnids is special as it accepts data values as well. A data value can be specified after each name.
Op	Highlight Operation
rnids	Highlight reaction names or frame ids.
rsubs	Highlight reaction substrings.
recns	Highlight reaction EC numbers.
pnids	Highlight pathway names or frame ids.
psubs	Highlight pathway substrings.
gnids	Highlight gene names or frame ids.
gsubs	Highlight gene substrings.
enids	Highlight enzyme names or frame ids.
esubs	Highlight enzyme substrings.
cnids	Highlight compound names or frame ids.
csubs	Highlight compound substrings.
xnids	Highlight a mix of names and frame ids with or without omics data.
pcids	Highlight pathways based on curation status.
rcids	Highlight reactions based on curation status.
revis	Highlight reactions based on evidence selected.
pevis	Highlight pathways based on evidence selected.
gregs	Highlight genes based on regulation selected.
greps	Highlight genes based on replicons.
The string specified after the ‘=’ for an operation must not be quoted, and any special characters must be URL encoded. The string is not case-sensitive. Multiple values separated by commas may be supplied for any of the above operations (except xnids, which uses URL-encoded spaces as separators). There is no difference between suppling multiple comma-separated values for a single operation versus specifying that operation multiple times, once for each value. When supplying omics data values using the xnids operation, the id and the value should be separated by a URL-encoded space. If multiple objects are specified using the xnids operation, either all must have data values (resulting in an omics viewer display), or any data values will be ignored (reverting to a basic highlighting operation).

Example highlighting URLs are as follows.


Examples of Cellular Overview Highlights
Example	URL
Highlight genes with a single color via op xnids	https://websvc.biocyc.org/overviewsWeb/celOv.shtml?zoomlevel=0&orgid=ECOLI&xnids=thiE trpA EG10790 b2518 b1704
Highlight genes with integer omics data via op xnids	https://websvc.biocyc.org/overviewsWeb/celOv.shtml?zoomlevel=0&orgid=ECOLI&xnids=thiE 4 trpA 6 EG10790 3 b2518 2 b1704 2.5
Highlight metabolites with float omics data via op xnids	https://websvc.biocyc.org/overviewsWeb/celOv.shtml?zoomlevel=1&orgid=ECOLI&xnids=L-BETA-ASPARTYL-P%204.0%20L-tryptophan%206.0%20L-LACTATE%203.0%20chorismate%202.0%20TYR%202.5 https://websvc.biocyc.org/overviewsWeb/celOv.shtml?zoomlevel=1&orgid=ECOLI&xnids=L-BETA-ASPARTYL-P 4.0 L-tryptophan 6.0 L-LACTATE 3.0 chorismate 2.0 TYR 2.5
Highlight reactions with omics data via op xnids	https://websvc.biocyc.org/overviewsWeb/celOv.shtml?zoomlevel=0&orgid=ECOLI&xnids=DIHYDROFOLATEREDUCT-RXN 4.0 DLACTDEHYDROGNAD-RXN 6.0 PGLUCISOM-RXN 3.0
Highlight reactions with omics data via op xnids	https://websvc.biocyc.org/overviewsWeb/celOv.shtml?zoomlevel=0&orgid=ECOLI&xnids=2.0 DIHYDROFOLATEREDUCT-RXN 4.0 DLACTDEHYDROGNAD-RXN&csubs=pyruv
Highlight enzyme based on substrings via op esubs	https://websvc.biocyc.org/overviewsWeb/celOv.shtml?zoomlevel=0&orgid=ECOLI&esubs=hydro&esubs=oxy&gsubs=arg,his
Highlight reactions based on evidence code via op revis	https://websvc.biocyc.org/overviewsWeb/celOv.shtml?zoomlevel=0&orgid=ECOLI&revis=EV-COMP
URLs with highlighting operations can be automatically generated by the Cellular Overview by applying the desired highlights using its interactive web interface, and then right-clicking on the diagram and selecting the command Generate Bookmark for Current Cellular Overview.

Paint Single Omics Data from a File on Metabolic Network Diagram
The metabolic network diagram can be invoked by using a URL link that specifies an omics data file that resides on an accessible website. See the Omics Viewer section of the Website Users Guide for more information about the Omics Viewer, including the data file format. This GET request will navigate to the Cellular Omics Viewer page and load it with data from the supplied URL and the parameters specified in the request itself.
URLs to paint omics data on the cellular overview diagram are of the form:

https://websvc.biocyc.org/overviewsWeb/celOv.shtml?orgid=[ORGID]&omics=t&zoomlevel=[0-3]&[PARAMETER]=[VALUE]

where

[ORGID] is the identifier for the organism database, e.g. ECOLI, BSUB, AFER243159
zoomlevel is an integer that indicates the extent to which the diagram has been zoomed. Zoomlevel 0 is the lowest magnification, level 3 is the largest magnification. If not supplied, the zoomlevel defaults to 0.
Multiple parameter/value pairs can be specified, but two are required:
url -- The URL of the omics data file
column1 -- An integer specifying the column number for the data from the omics data file that is to be displayed (column zero is the first column)
The list of available parameters are those marked O in the Omics Data Parameters Table.

Examples of Cellular Overview Single Omics Data Painting
Example	URL
Paint transcriptomics data from a sample data file. To supply your own data use a POST. Six columns of data are provided in this example; click triangle "Start" button to play the animation.	https://websvc.biocyc.org/overviewsWeb/celOv.shtml?orgid=ECOLI&omics=t&url=expr-examples/ecoli-feuer-toaerobic-significant.txt&expressiontype=relative&column1=1-6&log=on&class=gene
Paint metabolomics data from a sample data file. To supply your own data use a POST. Five columns of data are provided in this example; click triangle "Start" button to play the animation.	https://websvc.biocyc.org/overviewsWeb/celOv.shtml?orgid=HUMAN&omics=t&url=expr-examples/human-metabolomics-syn.txt&expressiontype=relative&column1=1-5&log=on&class=compound
Paint Multi Omics Data From File via POST
The multi omics viewer gives you the ability to upload up to four omics datasets onto the cellular overview. Each dataset is presented via a separate “visual channel.” The available channels are node (metabolite) colors, edge (reaction) colors, node size, and edge thickness. Typically, nodes are used to visualize metabolomics data and edges are used to visualize transcriptomics, proteomics, and reaction-flux data.

The cellular overview can be invoked and passed files via post by calling the overview-multi-omics-process-upr web service. One of the requirements of this post request is that the same session needs to be maintained.

Send a post request to: https://websvc.biocyc.org/overview-multi-omics-process-upr

POST Request Parameters:

Single File POST Request:
Should only be passed one file that would contain the master fields and omics data.
For more information relating to multi omics file formatting go here: link.

Parameters:

orgid:
Expects: The organism id
Required: yes

multifile:
To use single file mode must be false
Expects: t or f
Required: yes

file:
Expects: a file
Includes both the master fields and omics data.
Required: yes

Example:

With the parameters:

data = {
  "orgid" : "ECOLI",
  "multifile" : "f"
}
files = {  "file" : open("/file1/location/file1.txt "rb") }
Single file example file:
multi-omics-single-file-example.txt

Multi File POST Request:
Will be passed up to five files containing the master fields and omics data.

For more information relating to multi omics file formatting go here: link.

Parameters:

orgid:
Expects: The organism id
Required: yes

multifile:
To use multi file mode, must be true
Expects: t or f
Required: yes

masterfile:
Expects: A master file which includes all of the metafields.
Includes only the master file
Required: yes

file1 - file4:
Expects: a file
Includes only omics data
Required: yes – at minimum file1 and file2

Example:

With the parameters:

data = {
  "orgid" : "ECOLI",
  "multifile" : "t" 
}

files = {
  "masterfile" : open("/file1/location/masterfile.txt "rb")
  "file1" : open("/file1/location/file1.txt "rb")
  "file2" : open("/file2/location/file2.txt "rb")
  "file3" : open("/file3/location/file3.txt "rb")
  "file4" : open("/file4/location/file4.txt "rb")
}
Multi Omics file example files:
multi-omics-masterfile-example.txt
multi-omics-file1-example.txt
multi-omics-file2-example.txt

Python Library Requirements
The following libraries are needed in order to use the multiomics celov POST request:

requests
selenium
webdriver_manager
Download requirements.txt

To install the packages via the python package manager pip:

pip install -r /path/to/requirements.txt
Python Script Explanation
The script is first creating a session by logging into biocyc. Then after it is successfull, using the same session it makes a POST request to "/overview-multi-omics-process-upr". Afterwords it takes the returned url response, opens up a new firefox window using selenium with your current session, and it will redirect to the cellular overview page which will then load the posted data into the diagram.

What variables need modifed:

login_info: Swap the placeholder values to your biocyc credentials

data: Modify the data variable containing the two parameters to switch if it is single or muilti file mode and which organism id you are selecting.

files: Then change the files variable to reference each file you need for the POST request.

Script for POST request
import requests
from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.firefox.service import Service
from selenium.webdriver.firefox.options import Options
from selenium.webdriver.firefox.firefox_profile import FirefoxProfile
from webdriver_manager.firefox import GeckoDriverManager
	
base_url = "https://websvc.biocyc.org/"
login_url = base_url + "credentials/login/"
login_info = {
	"email" : "your_email",
	"password" : "your_password"
}
	
with requests.Session() as session:
	resp = session.post(login_url, data=login_info)
	if resp.status_code == 200:
		firefox_options = Options()
		firefox_options.add_argument("--start-maximized")
		url =  base_url + "/overview-multi-omics-process-upr"
		data = { "orgid" : "ECOLI", "multifile" : "t"}
		files = {  "masterfile" : open("/masterfile/location/masterfile.txt", "rb"), "file1" : open("/file1/location/file1.txt", "rb") , "file2" : open("/file1/location/file1.txt"", "rb") }
		headers = {'Cache-Control' : 'no-cache, no-store, must-revalidate'}
		response = session.post(url, headers=headers, data=data, files=files)
		print (response.status_code)
		print (response)
		if response.status_code == 200:
			print(response.text)
			omicsurl = "https://" + response.text
			driver = webdriver.Firefox(service=Service(GeckoDriverManager().install()), options=firefox_options)
			driver.get(base_url)
			for cookie in session.cookies:
				if cookie.name == "PTools-session" or cookie.name == "userIdentifier" or cookie.name == "recentOrgID0":
					item = {
						'name': cookie.name,
						'value': cookie.value,
						'domain': cookie.domain,
						'path': cookie.path,
						'expiry': cookie.expires,
						'secure': cookie.secure
					}
					driver.add_cookie(item)
			driver.refresh()
				
			print(response.text)
			driver.get(omicsurl)
			input()
			driver.quit()
	else:
		print (f"failed to render, status code: {response.status_code}")
Download celov_multiomics_post.py

SmartTables
SmartTable Data Retrieval
SmartTable data can be accessed from a web service in the following formats: JSON, XML, or tab delimited. This web service will request your account username and password.

The basic URL to retrieve a SmartTable is:

https://websvc.biocyc.org/st-get?id=[SMARTTABLE-ID]&format=[json|xml|tsv]

where

[SMARTTABLE-ID] is the identifier of the SmartTable, e.g. biocyc11-2762-3605383628, biocyc10-NIL-3608637901
[json|xml|tsv] is the requested format for the SmartTable data. json will return the data in JSON format, xml in XML format, and tsv in tab delimited format.
An example of the JSON structure can be found in the JSON Structure for SmartTables Retrieval and Creation.

Example URL:

https://websvc.biocyc.org/st-service-get?id=biocyc14-15682-3673724563&format=xml
SmartTable Creation
To create a SmartTable from data supplied via a file using the web service URL, the user must use a tool that can peform a "PUT" method. The result returned will contain the identifier of the created SmartTable. curl will be used for examples in this documentation section. When using curl, the terminal command may look like the following:

curl -u '[USER@DOMAIN.COM]' -X PUT -T [ST-PATH] https://websvc.biocyc.org/st-create?format=[json|xml|tsv]&orgid=[ORGID]&class=[CLASS]

where

[USER@DOMAIN.COM] is the user email address associated with their account
[ST-PATH] is the local path to the file containing SmartTable data, e.g. /tmp/test-st
[json|xml|tsv] is the data format of the SmartTable data stored at [st-path]
[ORGID] is the identifier for the organism database (tsv only)
[CLASS] is the class to which the objects in the first column belong, e.g. Genes, Compounds (tsv only)
Example terminal command:

curl -u 'USER@DOMAIN.COM' -X PUT -T /tmp/test-st.xml http://www.biocyc.org/st-create?format=xml
For tsv (tab-separated-values) input files, the orgid and type args are required; for json and xml input files, those parameters are encoded as part of the file format. For tsv input files, the first row of the file is assumed to contain column headers.
JSON Structure for SmartTable Data Retrieval and Creation
JSON file structure:

name (optional) - name field of SmartTable
description (optional) - description field of SmartTable
Contents of the short form, which is for creating a single column of values only:
type: an object class that values will be coerced into for the given database into the SmartTable
pgdb: an org id for the database in which SmartTable objects will reside and be coerced with the given type
values: an array of strings, which will be coerced into objects for the SmartTable in a single column
Contents of the long form, which is for creating multiple columns of values:
columns: a list (JSON array) of columns, each of which is:
name (optional): a string that will be the name of the column
type: an object class that values will be coerced into for the SmartTable
rows: a JSON array, each element is a dictionary mapping column ids to values. Each value corresponds to one cell in the SmartTable and contains:
value: a string, number, or {frameid: , pgdb: } which will be coerced into an object for the SmartTable
JSON examples:

// Short form
{"name": "sample short-form group",
 "description": "sample short-form description",
 "pgdb": "ECOLI",
 "type": "Genes",
 "values": ["trpA", "trpB"]
}

// long form (2 columns)
{"name": "sample long-form group",
 "description": "This is a longer form of JSON used with a SmartTable.",
 "columns": [{"name": "Gene",
                      "type": "Genes"},
                    {"name": "Expression Level"}]
 "rows": [{{"frameid": "EG11204", "pgdb": "ECOLI"},
                 "1.3"}
      {{"frameid": "EG11205", "pgdb": "ECOLI"},
        "-2.2"}]
}
XML Structure for SmartTable Data Retrieval and Creation
The XML file structure follows a similar structure to the JSON format. The following is the XML schema:

<xs:schema attributeFormDefault="unqualified" elementFormDefault="qualified" xmlns:xs="http://www.w3.org/2001/XMLSchema">
  <xs:element name="GROUP">
    <xs:complexType>
      <xs:attribute type="xs:string" name="ID"/>
      <xs:attribute type="xs:string" name="NAME"/>
      <xs:sequence>
        <xs:element name="COLUMNS">
          <xs:complexType>
            <xs:sequence>
              <xs:element name="COLUMN" maxOccurs="unbounded" minOccurs="1">
                <xs:complexType>
                  <xs:simpleContent>
                    <xs:extension base="xs:string">
                      <xs:attribute type="xs:string" name="NAME" use="optional"/>
                      <xs:attribute type="xs:string" name="TYPE" use="optional"/>
                    </xs:extension>
                  </xs:simpleContent>
                </xs:complexType>
              </xs:element>
            </xs:sequence>
          </xs:complexType>
        </xs:element>
        <xs:element name="ROWS">
          <xs:complexType>
            <xs:sequence>
              <xs:element name="ROW" maxOccurs="unbounded" minOccurs="0">
                <xs:complexType>
                  <xs:sequence>
                    <xs:element name="CELL" maxOccurs="unbounded" minOccurs="0">
                      <xs:complexType mixed="true">
                        <xs:sequence>
                          <xs:element name="FRAME" maxOccurs="unbounded" minOccurs="0">
                            <xs:complexType>
                              <xs:simpleContent>
                                <xs:extension base="xs:string">
                                  <xs:attribute type="xs:string" name="ID" use="optional"/>
                                  <xs:attribute type="xs:string" name="PGDB" use="optional"/>
                                </xs:extension>
                              </xs:simpleContent>
                            </xs:complexType>
                          </xs:element>
                          <xs:element type="xs:string" name="VALUE" maxOccurs="unbounded" minOccurs="0"/>
                        </xs:sequence>
                      </xs:complexType>
                    </xs:element>
                  </xs:sequence>
                </xs:complexType>
              </xs:element>
            </xs:sequence>
          </xs:complexType>
        </xs:element>
      </xs:sequence>
    </xs:complexType>
  </xs:element>
</xs:schema>
XML example:

<GROUP ID="biocyc14-2762-3579243971" NAME="Import from snp.txt">

  <COLUMNS>
    <COLUMN NAME="Position" TYPE="NUMBER"/>
    <COLUMN NAME="Gene" TYPE="All-Genes"/>

    <COLUMN
     NAME="Map to Alicyclobacillus acidocaldarius acidocaldarius DSM 446"
     TYPE="All-Genes"/>
    <COLUMN NAME="new column"/>

    <COLUMN NAME="Binding sites upstream of gene" TYPE="DNA-Binding-Sites"/>
  </COLUMNS>

  <ROWS>

    <ROW>
      <CELL><VALUE>10</VALUE></CELL>
      <CELL><FRAME ID="EG11024" PGDB="ECOLI"/></CELL>
      <CELL><FRAME ID="GCIO-1719" PGDB="AACI521098"/></CELL>
      <CELL/>

      <CELL>
        <FRAME ID="BS00206" PGDB="ECOLI"/>
        <FRAME ID="BS0-4061" PGDB="ECOLI"/>
        <FRAME ID="BS0-4062" PGDB="ECOLI"/>
      </CELL>
    </ROW>
  </ROWS>
</GROUP>
Other SmartTable Web Services
Various operations to manipulate existing SmartTables can be done with the following web services. Unless otherwise indicated, all use the HTTP GET method. In all example URLs, replace "[ID]" with the SmartTable identifier.

Service Name	Arguments	Example	Description
st-transform	id, transformid, index
https://websvc.biocyc.org/st-transform?id=[ID]&transformid=compound-producing-pathways&index=0	Adds a transformation column to an existing SmartTable.
id: SmartTable identifier.
transformid: Identifier of the transformation to use. See table of available transformations.
index: Column number to base the new transformation column on. First column is index 0.
st-property	id, propertyid, index
https://websvc.biocyc.org/st-property?id=[ID]&propertyid=common-name&index=0	Adds a slot property column to an existing SmartTable.
id: SmartTable identifier.
propertyid: Name of the slot to use.
index: Column number to base the new transformation column on. First column is index 0.
st-delete-rows	id, indices
https://websvc.biocyc.org/st-delete-rows?id=[ID]&indices=0,1	Deletes rows from an existing SmartTable.
id: SmartTable identifier.
indices: Comma-separated list of row numbers to delete. First row is index 0.
st-delete-column	id, index
https://websvc.biocyc.org/st-delete-column?id=[ID]&index=1	Deletes a column from an existing SmartTable.
id: SmartTable identifier
index: Column number to delete. First column is index 0.
st-modify-name	id, name
https://websvc.biocyc.org/st-modify-name?id=[ID]&name=All+Genes	Modifies the name of an existing SmartTable.
id: SmartTable identifier
name: SmartTable name text, url-encoded.
st-modify-description	id, description
https://websvc.biocyc.org/st-modify-description?id=[ID]&description=All+Genes+of+Dataset	Modifies the description of an existing SmartTable.
id: SmartTable identifier
description: SmartTable description text, url-encoded.
st-copy	id
https://websvc.biocyc.org/st-copy?id=[ID]	Creates a copy of an existing SmartTable.
id: SmartTable identifier
st-delete	id
https://websvc.biocyc.org/st-delete?id=[ID]	Deletes an existing SmartTable.
id: SmartTable identifier
st-enrichment	With GET: id, key, type, threshold, statistic, correction

With PUT: format, orgid, class, key, type, threshold, statistic, correction	https://websvc.biocyc.org/st-enrichment?id=[ID]&type=enrichment&key=enrich-genes-pwys	With the GET method, runs an enrichment analysis on the first column of an existing SmartTable. To upload a file of genes or compounds and run an enrichment analysis on it in one step, use the PUT method, and also supply the args needed for the st-create service.
id: SmartTable identifier (GET only)
key: one value from enrich-all-gen, enrich-go-cc, enrich-go-bp, enrich-go-mf, enrich-genes-regulators-indirect, enrich-genes-regulators, enrich-cpds-pwys, or enrich-genes-pwys
type: enrichment, depletion, or enrichment-and-depletion
threshold: a number between 0 and 1; only include results with P-value less than this number (default = 0.1)
statistic: fisher-exact (default), parent-child-union, or parent-child-intersection
correction: none (default), bonferroni, bh (Benjamini-Hochberg), or by (Benjamini-Yekutieli)
format: json, xml or tsv; the format of the uploaded data (PUT only)
class: Genes or Compounds (PUT, format=tsv only)
orgid: the identifier for the organism database (PUT, format=tsv only)