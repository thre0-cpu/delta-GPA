BRENDA Enzyme Database
SOAP access
Use of this online version of BRENDA is free under the CC BY 4.0 license. See terms of use for full details. Please register to use BRENDA's web service!

In the following, the source code of complete SOAP clients are listed for the programming languages Perl, PHP, Python and Java.
They represent example clients for extracting all KM_Value entities of BRENDA - in this example representing the EC number '1.1.1.1' and the organism 'Homo sapiens' - by using the corresponding method getKmValue(String).
In order to adapt these SOAP clients for other SOAP methods, only the marked yellow lines of source code have to be replaced by the code snippets listed under the respective method (see below). In order to use BRENDA's web service you need a valid email address and password.
Content
SOAP client examples

SOAP methods
Ligand structure Id
Reference by Id
Reference by pubmed Id
Activating Compound
Application
CAS Registry Number
Cloned
Cofactor
Crystallization
Disease
EC Number
Protein Variants (Engineering)
Enzyme Names
 
Expression
General Information
General Stability
IC50 Value
Inhibitors
KCat KM Value
KI Value
KM Value
Ligands
Localization
Metals Ions
Molecular Weight
Natural Product
 
Natural Substrate
Natural Substrates Products
Organic Solvent Stability
Organism
Oxidation Stability
Pathway
PDB
pH Optimum
pH Range
pH Stability
pI Value
Posttranslational Modification
Product
 
Purification
Reaction
Reaction Type
Recommended Name
Reference
Renatured
Sequence
Source Tissue
Specific Activity
Storage Stability
Substrate
Substrates Products
Subunits
 
Synonyms
Systematic Name
Temperature Optimum
Temperature Range
Temperature Stability
Turnover Number

SOAP client examples
Python 3
prerequisite: Installation of Zeep
#!/usr/bin/python
from zeep import Client, Settings
import hashlib

wsdl = "https://www.brenda-enzymes.org/soap/brenda_zeep.wsdl"
password = hashlib.sha256("myPassword".encode("utf-8")).hexdigest()
settings = Settings(strict=False)
client = Client(wsdl, settings=settings)
parameters = ( "j.doe@example.edu",password,"ecNumber*1.1.1.1","organism*Homo sapiens","kmValue*",
              "kmValueMaximum*","substrate*","commentary*","ligandStructureId*","literature*" )
resultString = client.service.getKmValue(*parameters)
print (resultString)

SOAP Methods
Ligand structure ID
1. getLigandStructureIdByCompoundName(string)
Input
top
Python3 example code snippet:
parameters = ("j.doe@example.edu",password,"Zn2+")
resultString = client.service.getLigandStructureIdByCompoundName(*parameters)
Output
String containing the ligand structure Id (BRENDA group ID), e.g. "28905"

Reference by Id
2. getReferenceById(String)
Input
top
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password,"64375"
resultString = client.service.getReferenceById(*parameters)
Output
String containing all fields of the specified Reference (fields separated by #) or a reference object (Python 3), e.g."authors*string#title*string#journal*string#volume*string#pages*string#year*string#pubmedId*string#" or "{'authors':string, 'title':string, 'journal':string, 'volume'"string, 'pages':string, 'year':string, 'pubmedId':string}"

Reference by Pubmed Id
3. getReferenceByPubmedId(String)
Input
top
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password,"9514116"
resultString = client.service.getReferenceByPubmedId(*parameters)
Output
String containing all fields of the specified Reference (fields separated by #) or a reference object, e.g. "authors*string#title*string#journal*string#volume*string#pages*string#year*string#pubmedId*string#" or "{'authors':string, 'title':string, 'journal':string, 'volume'"string, 'pages':string, 'year':string, 'pubmedId':string}"

Please consider: About 12% of the more than 100,000 references contained in BRENDA do not possess a PubMed ID because the corresponding articles/journals are not of biomedical relevance. Therefore, they are not listed in PubMed and, thus, these articles are not accessible via a PubMed ID. Please use either the method getReferenceById(String referenceID) (using the BRENDA reference ID) or the method getReference(String) instead.
Activating Compound
4. getEcNumbersFromActivatingCompound()
Input
top

Input Data Type String (for all languages):
email,password
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password)
resultString = client.service.getEcNumbersFromActivatingCompound(*parameters)
Output
String containing the different EC Numbers (separated by !) or an array of the different EC Numbers (Python 3) linked to Activating Compound, e.g. "EC Number1!EC Number2!EC Number3!EC Number4!..." or "['EC Number1','EC Number2','EC Number3','EC Number4']"

5. getOrganismsFromActivatingCompound()
Input
top
Input Data Type String (for all languages):
email,password
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password,"ecNumber*", "activatingCompound*", "commentary*", "organism*", "ligandStructureId*", "literature*")
resultString = client.service.getOrganismsFromActivatingCompound(*parameters)
Output
String containing the different Organisms (separated by !) or an array of Organisms (Python 3) linked to Activating Compound, e.g. "Organism1!Organism2!Organism3!Organism4!..." or "['Organism1','Organism2','Organism3','Organism4']"

6. getActivatingCompound(String)
Input
top
Input Data Type String (for all languages):

email,password,"ecNumber*1.1.1.1#organism*Mus musculus"
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password,"ecNumber*", "activatingCompound*", "commentary*", "organism*", "ligandStructureId*", "literature*")
resultString = client.service.getActivatingCompound(*parameters)
Either the key field ecNumber (e.g. "ecNumber*1.1.1.2"), the key field organism (e.g. "organism*Homo sapiens") or both (e.g. "ecNumber*1.1.1.2#organism*Homo sapiens") have to be specified. If none of these key fields is used as an input parameter, an empty String will be returned by the SOAP query.
In addition, the following optional parameters can be specified: "activatingCompound", "commentary", "ligandStructureId", "literature
Output
String containing the Activating Compound entries (entries separated by !, fields separated by #) or an array of Activating Compound objects (Python 3), e.g. "ecNumber*string#activatingCompound*string#commentary*string#organism*string#ligandStructureId*string!
 ecNumber*string#activatingCompound*string#commentary*string#organism*string#ligandStructureId*string!..."
or
[{'ecNumber':string, 'activatingCompound':string, 'commentary':string, 'organism':string, 'ligandStructureId':string, 'literature':[int1,int2,int3,...]},{'ecNumber':string, 'activatingCompound':string, 'commentary':string, 'organism':string, 'ligandStructureId':string, 'literature':[int1,int2,int3,...]},...]"

Application
7. getEcNumbersFromApplication()
Input
top

Input Data Type String (for all languages):
email,password
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password)
resultString = client.service.getEcNumbersFromApplication(*parameters)
Output
String containing the different EC Numbers (separated by !) or an array of the different EC Numbers (Python 3) linked to Application, e.g. "EC Number1!EC Number2!EC Number3!EC Number4!..." or "['EC Number1','EC Number2','EC Number3','EC Number4']"

8. getOrganismsFromApplication()
Input
top
Input Data Type String (for all languages):
email,password
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password,"ecNumber*", "application*", "commentary*", "organism*", "literature*", "ecNumber*")
resultString = client.service.getOrganismsFromApplication(*parameters)
Output
String containing the different Organisms (separated by !) or an array of Organisms (Python 3) linked to Application, e.g. "Organism1!Organism2!Organism3!Organism4!..." or "['Organism1','Organism2','Organism3','Organism4']"

9. getApplication(String)
Input
top
Input Data Type String (for all languages):

email,password,"ecNumber*1.1.1.1#organism*Mus musculus"
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password,"ecNumber*", "application*", "commentary*", "organism*", "literature*", "ecNumber*")
resultString = client.service.getApplication(*parameters)
Either the key field ecNumber (e.g. "ecNumber*1.1.1.2"), the key field organism (e.g. "organism*Homo sapiens") or both (e.g. "ecNumber*1.1.1.2#organism*Homo sapiens") have to be specified. If none of these key fields is used as an input parameter, an empty String will be returned by the SOAP query.
In addition, the following optional parameters can be specified: "application", "commentary", "literature
Output
String containing the Application entries (entries separated by !, fields separated by #) or an array of Application objects (Python 3), e.g. "ecNumber*string#application*string#commentary*string#organism*string#ecNumber*string!
 ecNumber*string#application*string#commentary*string#organism*string#ecNumber*string!..."
or
[{'ecNumber':string, 'application':string, 'commentary':string, 'organism':string, 'literature':[int1,int2,int3,...], 'ecNumber':string},{'ecNumber':string, 'application':string, 'commentary':string, 'organism':string, 'literature':[int1,int2,int3,...], 'ecNumber':string},...]"

CAS Registry Number
10. getEcNumbersFromCasRegistryNumber()
Input
top

Input Data Type String (for all languages):
email,password
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password)
resultString = client.service.getEcNumbersFromCasRegistryNumber(*parameters)
Output
String containing the different EC Numbers (separated by !) or an array of the different EC Numbers (Python 3) linked to CAS Registry Number, e.g. "EC Number1!EC Number2!EC Number3!EC Number4!..." or "['EC Number1','EC Number2','EC Number3','EC Number4']"

11. getCasRegistryNumber(String)
Input
top
Input Data Type String (for all languages):

email,password,"ecNumber*1.1.1.1"
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password,"ecNumber*", "casRegistryNumber*", "commentary*")
resultString = client.service.getCasRegistryNumber(*parameters)
Either the key field ecNumber (e.g. "ecNumber*1.1.1.2"), the key field organism (e.g. "organism*Homo sapiens") or both (e.g. "ecNumber*1.1.1.2#organism*Homo sapiens") have to be specified. If none of these key fields is used as an input parameter, an empty String will be returned by the SOAP query.
In addition, the following optional parameters can be specified: "casRegistryNumber", "commentary
Output
String containing the CAS Registry Number entries (entries separated by !, fields separated by #) or an array of CAS Registry Number objects (Python 3), e.g. "ecNumber*string#casRegistryNumber*string#commentary*string!
 ecNumber*string#casRegistryNumber*string#commentary*string!..."
or
[{'ecNumber':string, 'casRegistryNumber':string, 'commentary':string},{'ecNumber':string, 'casRegistryNumber':string, 'commentary':string},...]"

Cloned
12. getEcNumbersFromCloned()
Input
top

Input Data Type String (for all languages):
email,password
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password)
resultString = client.service.getEcNumbersFromCloned(*parameters)
Output
String containing the different EC Numbers (separated by !) or an array of the different EC Numbers (Python 3) linked to Cloned, e.g. "EC Number1!EC Number2!EC Number3!EC Number4!..." or "['EC Number1','EC Number2','EC Number3','EC Number4']"

13. getOrganismsFromCloned()
Input
top
Input Data Type String (for all languages):
email,password
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password,"ecNumber*", "commentary*", "organism*", "literature*")
resultString = client.service.getOrganismsFromCloned(*parameters)
Output
String containing the different Organisms (separated by !) or an array of Organisms (Python 3) linked to Cloned, e.g. "Organism1!Organism2!Organism3!Organism4!..." or "['Organism1','Organism2','Organism3','Organism4']"

14. getCloned(String)
Input
top
Input Data Type String (for all languages):

email,password,"ecNumber*1.1.1.1#organism*Mus musculus"
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password,"ecNumber*", "commentary*", "organism*", "literature*")
resultString = client.service.getCloned(*parameters)
Either the key field ecNumber (e.g. "ecNumber*1.1.1.2"), the key field organism (e.g. "organism*Homo sapiens") or both (e.g. "ecNumber*1.1.1.2#organism*Homo sapiens") have to be specified. If none of these key fields is used as an input parameter, an empty String will be returned by the SOAP query.
In addition, the following optional parameters can be specified: "commentary", "literature
Output
String containing the Cloned entries (entries separated by !, fields separated by #) or an array of Cloned objects (Python 3), e.g. "ecNumber*string#commentary*string#organism*string!
 ecNumber*string#commentary*string#organism*string!..."
or
[{'ecNumber':string, 'commentary':string, 'organism':string, 'literature':[int1,int2,int3,...]},{'ecNumber':string, 'commentary':string, 'organism':string, 'literature':[int1,int2,int3,...]},...]"

Cofactor
15. getEcNumbersFromCofactor()
Input
top

Input Data Type String (for all languages):
email,password
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password)
resultString = client.service.getEcNumbersFromCofactor(*parameters)
Output
String containing the different EC Numbers (separated by !) or an array of the different EC Numbers (Python 3) linked to Cofactor, e.g. "EC Number1!EC Number2!EC Number3!EC Number4!..." or "['EC Number1','EC Number2','EC Number3','EC Number4']"

16. getOrganismsFromCofactor()
Input
top
Input Data Type String (for all languages):
email,password
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password,"ecNumber*", "cofactor*", "commentary*", "organism*", "ligandStructureId*", "literature*")
resultString = client.service.getOrganismsFromCofactor(*parameters)
Output
String containing the different Organisms (separated by !) or an array of Organisms (Python 3) linked to Cofactor, e.g. "Organism1!Organism2!Organism3!Organism4!..." or "['Organism1','Organism2','Organism3','Organism4']"

17. getCofactor(String)
Input
top
Input Data Type String (for all languages):

email,password,"ecNumber*1.1.1.1#organism*Mus musculus"
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password,"ecNumber*", "cofactor*", "commentary*", "organism*", "ligandStructureId*", "literature*")
resultString = client.service.getCofactor(*parameters)
Either the key field ecNumber (e.g. "ecNumber*1.1.1.2"), the key field organism (e.g. "organism*Homo sapiens") or both (e.g. "ecNumber*1.1.1.2#organism*Homo sapiens") have to be specified. If none of these key fields is used as an input parameter, an empty String will be returned by the SOAP query.
In addition, the following optional parameters can be specified: "cofactor", "commentary", "ligandStructureId", "literature
Output
String containing the Cofactor entries (entries separated by !, fields separated by #) or an array of Cofactor objects (Python 3), e.g. "ecNumber*string#cofactor*string#commentary*string#organism*string#ligandStructureId*string!
 ecNumber*string#cofactor*string#commentary*string#organism*string#ligandStructureId*string!..."
or
[{'ecNumber':string, 'cofactor':string, 'commentary':string, 'organism':string, 'ligandStructureId':string, 'literature':[int1,int2,int3,...]},{'ecNumber':string, 'cofactor':string, 'commentary':string, 'organism':string, 'ligandStructureId':string, 'literature':[int1,int2,int3,...]},...]"

Crystallization
18. getEcNumbersFromCrystallization()
Input
top

Input Data Type String (for all languages):
email,password
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password)
resultString = client.service.getEcNumbersFromCrystallization(*parameters)
Output
String containing the different EC Numbers (separated by !) or an array of the different EC Numbers (Python 3) linked to Crystallization, e.g. "EC Number1!EC Number2!EC Number3!EC Number4!..." or "['EC Number1','EC Number2','EC Number3','EC Number4']"

19. getOrganismsFromCrystallization()
Input
top
Input Data Type String (for all languages):
email,password
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password,"ecNumber*", "commentary*", "organism*", "literature*")
resultString = client.service.getOrganismsFromCrystallization(*parameters)
Output
String containing the different Organisms (separated by !) or an array of Organisms (Python 3) linked to Crystallization, e.g. "Organism1!Organism2!Organism3!Organism4!..." or "['Organism1','Organism2','Organism3','Organism4']"

20. getCrystallization(String)
Input
top
Input Data Type String (for all languages):

email,password,"ecNumber*1.1.1.1#organism*Mus musculus"
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password,"ecNumber*", "commentary*", "organism*", "literature*")
resultString = client.service.getCrystallization(*parameters)
Either the key field ecNumber (e.g. "ecNumber*1.1.1.2"), the key field organism (e.g. "organism*Homo sapiens") or both (e.g. "ecNumber*1.1.1.2#organism*Homo sapiens") have to be specified. If none of these key fields is used as an input parameter, an empty String will be returned by the SOAP query.
In addition, the following optional parameters can be specified: "commentary", "literature
Output
String containing the Crystallization entries (entries separated by !, fields separated by #) or an array of Crystallization objects (Python 3), e.g. "ecNumber*string#commentary*string#organism*string!
 ecNumber*string#commentary*string#organism*string!..."
or
[{'ecNumber':string, 'commentary':string, 'organism':string, 'literature':[int1,int2,int3,...]},{'ecNumber':string, 'commentary':string, 'organism':string, 'literature':[int1,int2,int3,...]},...]"

Disease
21. getEcNumbersFromDisease()
Input
top

Input Data Type String (for all languages):
email,password
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password)
resultString = client.service.getEcNumbersFromDisease(*parameters)
Output
String containing the different EC Numbers (separated by !) or an array of the different EC Numbers (Python 3) linked to Disease, e.g. "EC Number1!EC Number2!EC Number3!EC Number4!..." or "['EC Number1','EC Number2','EC Number3','EC Number4']"

22. getDisease(String)
Input
top
Input Data Type String (for all languages):

email,password,"ecNumber*1.1.1.1"
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password,"ecNumber*", "disease*", "pubmedId*", "titlePub*", "category*", "highestConfidenceLevel*")
resultString = client.service.getDisease(*parameters)
Either the key field ecNumber (e.g. "ecNumber*1.1.1.2"), the key field organism (e.g. "organism*Homo sapiens") or both (e.g. "ecNumber*1.1.1.2#organism*Homo sapiens") have to be specified. If none of these key fields is used as an input parameter, an empty String will be returned by the SOAP query.
In addition, the following optional parameters can be specified: "disease", "pubmedId", "titlePub", "category", "highestConfidenceLevel
Output
String containing the Disease entries (entries separated by !, fields separated by #) or an array of Disease objects (Python 3), e.g. "ecNumber*string#disease*string#pubmedId*string#titlePub*string#category*string#highestConfidenceLevel*string!
 ecNumber*string#disease*string#pubmedId*string#titlePub*string#category*string#highestConfidenceLevel*string!..."
or
[{'ecNumber':string, 'disease':string, 'pubmedId':string, 'titlePub':string, 'category':string, 'highestConfidenceLevel':string},{'ecNumber':string, 'disease':string, 'pubmedId':string, 'titlePub':string, 'category':string, 'highestConfidenceLevel':string},...]"

EC Number
23. getEcNumbersFromEcNumber()
Input
top

Input Data Type String (for all languages):
email,password
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password)
resultString = client.service.getEcNumbersFromEcNumber(*parameters)
Output
String containing the different EC Numbers (separated by !) or an array of the different EC Numbers (Python 3) linked to EC Number, e.g. "EC Number1!EC Number2!EC Number3!EC Number4!..." or "['EC Number1','EC Number2','EC Number3','EC Number4']"

24. getEcNumber(String)
Input
top
Input Data Type String (for all languages):

email,password,"ecNumber*1.1.1.1"
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password,"ecNumber*", "commentary*")
resultString = client.service.getEcNumber(*parameters)
Either the key field ecNumber (e.g. "ecNumber*1.1.1.2"), the key field organism (e.g. "organism*Homo sapiens") or both (e.g. "ecNumber*1.1.1.2#organism*Homo sapiens") have to be specified. If none of these key fields is used as an input parameter, an empty String will be returned by the SOAP query.
In addition, the following optional parameters can be specified: "commentary
Output
String containing the EC Number entries (entries separated by !, fields separated by #) or an array of EC Number objects (Python 3), e.g. "ecNumber*string#commentary*string!
 ecNumber*string#commentary*string!..."
or
[{'ecNumber':string, 'commentary':string},{'ecNumber':string, 'commentary':string},...]"

Protein Variants (Engineering)
25. getEcNumbersFromEngineering()
Input
top

Input Data Type String (for all languages):
email,password
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password)
resultString = client.service.getEcNumbersFromEngineering(*parameters)
Output
String containing the different EC Numbers (separated by !) or an array of the different EC Numbers (Python 3) linked to Engineering, e.g. "EC Number1!EC Number2!EC Number3!EC Number4!..." or "['EC Number1','EC Number2','EC Number3','EC Number4']"

26. getOrganismsFromEngineering()
Input
top
Input Data Type String (for all languages):
email,password
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password,"ecNumber*", "engineering*", "commentary*", "organism*", "literature*")
resultString = client.service.getOrganismsFromEngineering(*parameters)
Output
String containing the different Organisms (separated by !) or an array of Organisms (Python 3) linked to Engineering, e.g. "Organism1!Organism2!Organism3!Organism4!..." or "['Organism1','Organism2','Organism3','Organism4']"

27. getEngineering(String)
Input
top
Input Data Type String (for all languages):

email,password,"ecNumber*1.1.1.1#organism*Mus musculus"
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password,"ecNumber*", "engineering*", "commentary*", "organism*", "literature*")
resultString = client.service.getEngineering(*parameters)
Either the key field ecNumber (e.g. "ecNumber*1.1.1.2"), the key field organism (e.g. "organism*Homo sapiens") or both (e.g. "ecNumber*1.1.1.2#organism*Homo sapiens") have to be specified. If none of these key fields is used as an input parameter, an empty String will be returned by the SOAP query.
In addition, the following optional parameters can be specified: "engineering", "commentary", "literature
Output
String containing the Engineering entries (entries separated by !, fields separated by #) or an array of Engineering objects (Python 3), e.g. "ecNumber*string#engineering*string#commentary*string#organism*string!
 ecNumber*string#engineering*string#commentary*string#organism*string!..."
or
[{'ecNumber':string, 'engineering':string, 'commentary':string, 'organism':string, 'literature':[int1,int2,int3,...]},{'ecNumber':string, 'engineering':string, 'commentary':string, 'organism':string, 'literature':[int1,int2,int3,...]},...]"

Enzyme Names
28. getEcNumbersFromEnzymeNames()
Input
top

Input Data Type String (for all languages):
email,password
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password)
resultString = client.service.getEcNumbersFromEnzymeNames(*parameters)
Output
String containing the different EC Numbers (separated by !) or an array of the different EC Numbers (Python 3) linked to Enzyme Names, e.g. "EC Number1!EC Number2!EC Number3!EC Number4!..." or "['EC Number1','EC Number2','EC Number3','EC Number4']"

29. getEnzymeNames(String)
Input
top
Input Data Type String (for all languages):

email,password,"ecNumber*1.1.1.1"
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password,"ecNumber*", "synonyms*")
resultString = client.service.getEnzymeNames(*parameters)
Either the key field ecNumber (e.g. "ecNumber*1.1.1.2"), the key field organism (e.g. "organism*Homo sapiens") or both (e.g. "ecNumber*1.1.1.2#organism*Homo sapiens") have to be specified. If none of these key fields is used as an input parameter, an empty String will be returned by the SOAP query.
In addition, the following optional parameters can be specified: "synonyms
Output
String containing the Enzyme Names entries (entries separated by !, fields separated by #) or an array of Enzyme Names objects (Python 3), e.g. "ecNumber*string#synonyms*string!
 ecNumber*string#synonyms*string!..."
or
[{'ecNumber':string, 'synonyms':string},{'ecNumber':string, 'synonyms':string},...]"

Expression
30. getEcNumbersFromExpression()
Input
top

Input Data Type String (for all languages):
email,password
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password)
resultString = client.service.getEcNumbersFromExpression(*parameters)
Output
String containing the different EC Numbers (separated by !) or an array of the different EC Numbers (Python 3) linked to Expression, e.g. "EC Number1!EC Number2!EC Number3!EC Number4!..." or "['EC Number1','EC Number2','EC Number3','EC Number4']"

31. getOrganismsFromExpression()
Input
top
Input Data Type String (for all languages):
email,password
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password,"ecNumber*", "expression*", "commentary*", "organism*", "literature*", "ecNumber*")
resultString = client.service.getOrganismsFromExpression(*parameters)
Output
String containing the different Organisms (separated by !) or an array of Organisms (Python 3) linked to Expression, e.g. "Organism1!Organism2!Organism3!Organism4!..." or "['Organism1','Organism2','Organism3','Organism4']"

32. getExpression(String)
Input
top
Input Data Type String (for all languages):

email,password,"ecNumber*1.1.1.1#organism*Mus musculus"
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password,"ecNumber*", "expression*", "commentary*", "organism*", "literature*", "ecNumber*")
resultString = client.service.getExpression(*parameters)
Either the key field ecNumber (e.g. "ecNumber*1.1.1.2"), the key field organism (e.g. "organism*Homo sapiens") or both (e.g. "ecNumber*1.1.1.2#organism*Homo sapiens") have to be specified. If none of these key fields is used as an input parameter, an empty String will be returned by the SOAP query.
In addition, the following optional parameters can be specified: "expression", "commentary", "literature
Output
String containing the Expression entries (entries separated by !, fields separated by #) or an array of Expression objects (Python 3), e.g. "ecNumber*string#expression*string#commentary*string#organism*string#ecNumber*string!
 ecNumber*string#expression*string#commentary*string#organism*string#ecNumber*string!..."
or
[{'ecNumber':string, 'expression':string, 'commentary':string, 'organism':string, 'literature':[int1,int2,int3,...], 'ecNumber':string},{'ecNumber':string, 'expression':string, 'commentary':string, 'organism':string, 'literature':[int1,int2,int3,...], 'ecNumber':string},...]"

General Information
33. getEcNumbersFromGeneralInformation()
Input
top

Input Data Type String (for all languages):
email,password
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password)
resultString = client.service.getEcNumbersFromGeneralInformation(*parameters)
Output
String containing the different EC Numbers (separated by !) or an array of the different EC Numbers (Python 3) linked to General Information, e.g. "EC Number1!EC Number2!EC Number3!EC Number4!..." or "['EC Number1','EC Number2','EC Number3','EC Number4']"

34. getOrganismsFromGeneralInformation()
Input
top
Input Data Type String (for all languages):
email,password
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password,"ecNumber*", "generalInformation*", "commentary*", "organism*", "literature*")
resultString = client.service.getOrganismsFromGeneralInformation(*parameters)
Output
String containing the different Organisms (separated by !) or an array of Organisms (Python 3) linked to General Information, e.g. "Organism1!Organism2!Organism3!Organism4!..." or "['Organism1','Organism2','Organism3','Organism4']"

35. getGeneralInformation(String)
Input
top
Input Data Type String (for all languages):

email,password,"ecNumber*1.1.1.1#organism*Mus musculus"
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password,"ecNumber*", "generalInformation*", "commentary*", "organism*", "literature*")
resultString = client.service.getGeneralInformation(*parameters)
Either the key field ecNumber (e.g. "ecNumber*1.1.1.2"), the key field organism (e.g. "organism*Homo sapiens") or both (e.g. "ecNumber*1.1.1.2#organism*Homo sapiens") have to be specified. If none of these key fields is used as an input parameter, an empty String will be returned by the SOAP query.
In addition, the following optional parameters can be specified: "generalInformation", "commentary", "literature
Output
String containing the General Information entries (entries separated by !, fields separated by #) or an array of General Information objects (Python 3), e.g. "ecNumber*string#generalInformation*string#commentary*string#organism*string!
 ecNumber*string#generalInformation*string#commentary*string#organism*string!..."
or
[{'ecNumber':string, 'generalInformation':string, 'commentary':string, 'organism':string, 'literature':[int1,int2,int3,...]},{'ecNumber':string, 'generalInformation':string, 'commentary':string, 'organism':string, 'literature':[int1,int2,int3,...]},...]"

General Stability
36. getEcNumbersFromGeneralStability()
Input
top

Input Data Type String (for all languages):
email,password
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password)
resultString = client.service.getEcNumbersFromGeneralStability(*parameters)
Output
String containing the different EC Numbers (separated by !) or an array of the different EC Numbers (Python 3) linked to General Stability, e.g. "EC Number1!EC Number2!EC Number3!EC Number4!..." or "['EC Number1','EC Number2','EC Number3','EC Number4']"

37. getOrganismsFromGeneralStability()
Input
top
Input Data Type String (for all languages):
email,password
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password,"ecNumber*", "generalStability*", "organism*", "literature*")
resultString = client.service.getOrganismsFromGeneralStability(*parameters)
Output
String containing the different Organisms (separated by !) or an array of Organisms (Python 3) linked to General Stability, e.g. "Organism1!Organism2!Organism3!Organism4!..." or "['Organism1','Organism2','Organism3','Organism4']"

38. getGeneralStability(String)
Input
top
Input Data Type String (for all languages):

email,password,"ecNumber*1.1.1.1#organism*Mus musculus"
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password,"ecNumber*", "generalStability*", "organism*", "literature*")
resultString = client.service.getGeneralStability(*parameters)
Either the key field ecNumber (e.g. "ecNumber*1.1.1.2"), the key field organism (e.g. "organism*Homo sapiens") or both (e.g. "ecNumber*1.1.1.2#organism*Homo sapiens") have to be specified. If none of these key fields is used as an input parameter, an empty String will be returned by the SOAP query.
In addition, the following optional parameters can be specified: "generalStability", "literature
Output
String containing the General Stability entries (entries separated by !, fields separated by #) or an array of General Stability objects (Python 3), e.g. "ecNumber*string#generalStability*string#organism*string!
 ecNumber*string#generalStability*string#organism*string!..."
or
[{'ecNumber':string, 'generalStability':string, 'organism':string, 'literature':[int1,int2,int3,...]},{'ecNumber':string, 'generalStability':string, 'organism':string, 'literature':[int1,int2,int3,...]},...]"

IC50 Value
39. getEcNumbersFromIc50Value()
Input
top

Input Data Type String (for all languages):
email,password
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password)
resultString = client.service.getEcNumbersFromIc50Value(*parameters)
Output
String containing the different EC Numbers (separated by !) or an array of the different EC Numbers (Python 3) linked to IC50 Value, e.g. "EC Number1!EC Number2!EC Number3!EC Number4!..." or "['EC Number1','EC Number2','EC Number3','EC Number4']"

40. getOrganismsFromIc50Value()
Input
top
Input Data Type String (for all languages):
email,password
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password,"ecNumber*", "ic50Value*", "ic50ValueMaximum*", "inhibitor*", "commentary*", "organism*", "ligandStructureId*", "literature*")
resultString = client.service.getOrganismsFromIc50Value(*parameters)
Output
String containing the different Organisms (separated by !) or an array of Organisms (Python 3) linked to IC50 Value, e.g. "Organism1!Organism2!Organism3!Organism4!..." or "['Organism1','Organism2','Organism3','Organism4']"

41. getIc50Value(String)
Input
top
Input Data Type String (for all languages):

email,password,"ecNumber*1.1.1.1#organism*Mus musculus"
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password,"ecNumber*", "ic50Value*", "ic50ValueMaximum*", "inhibitor*", "commentary*", "organism*", "ligandStructureId*", "literature*")
resultString = client.service.getIc50Value(*parameters)
Either the key field ecNumber (e.g. "ecNumber*1.1.1.2"), the key field organism (e.g. "organism*Homo sapiens") or both (e.g. "ecNumber*1.1.1.2#organism*Homo sapiens") have to be specified. If none of these key fields is used as an input parameter, an empty String will be returned by the SOAP query.
In addition, the following optional parameters can be specified: "ic50Value", "ic50ValueMaximum", "inhibitor", "commentary", "ligandStructureId", "literature
Output
String containing the IC50 Value entries (entries separated by !, fields separated by #) or an array of IC50 Value objects (Python 3), e.g. "ecNumber*string#ic50Value*string#ic50ValueMaximum*string#inhibitor*string#commentary*string#organism*string#ligandStructureId*string!
 ecNumber*string#ic50Value*string#ic50ValueMaximum*string#inhibitor*string#commentary*string#organism*string#ligandStructureId*string!..."
or
[{'ecNumber':string, 'ic50Value':string, 'ic50ValueMaximum':string, 'inhibitor':string, 'commentary':string, 'organism':string, 'ligandStructureId':string, 'literature':[int1,int2,int3,...]},{'ecNumber':string, 'ic50Value':string, 'ic50ValueMaximum':string, 'inhibitor':string, 'commentary':string, 'organism':string, 'ligandStructureId':string, 'literature':[int1,int2,int3,...]},...]"

Inhibitors
42. getEcNumbersFromInhibitors()
Input
top

Input Data Type String (for all languages):
email,password
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password)
resultString = client.service.getEcNumbersFromInhibitors(*parameters)
Output
String containing the different EC Numbers (separated by !) or an array of the different EC Numbers (Python 3) linked to Inhibitors, e.g. "EC Number1!EC Number2!EC Number3!EC Number4!..." or "['EC Number1','EC Number2','EC Number3','EC Number4']"

43. getOrganismsFromInhibitors()
Input
top
Input Data Type String (for all languages):
email,password
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password,"ecNumber*", "inhibitor*", "commentary*", "organism*", "ligandStructureId*", "literature*")
resultString = client.service.getOrganismsFromInhibitors(*parameters)
Output
String containing the different Organisms (separated by !) or an array of Organisms (Python 3) linked to Inhibitors, e.g. "Organism1!Organism2!Organism3!Organism4!..." or "['Organism1','Organism2','Organism3','Organism4']"

44. getInhibitors(String)
Input
top
Input Data Type String (for all languages):

email,password,"ecNumber*1.1.1.1#organism*Mus musculus"
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password,"ecNumber*", "inhibitor*", "commentary*", "organism*", "ligandStructureId*", "literature*")
resultString = client.service.getInhibitors(*parameters)
Either the key field ecNumber (e.g. "ecNumber*1.1.1.2"), the key field organism (e.g. "organism*Homo sapiens") or both (e.g. "ecNumber*1.1.1.2#organism*Homo sapiens") have to be specified. If none of these key fields is used as an input parameter, an empty String will be returned by the SOAP query.
In addition, the following optional parameters can be specified: "inhibitor", "commentary", "ligandStructureId", "literature
Output
String containing the Inhibitors entries (entries separated by !, fields separated by #) or an array of Inhibitors objects (Python 3), e.g. "ecNumber*string#inhibitor*string#commentary*string#organism*string#ligandStructureId*string!
 ecNumber*string#inhibitor*string#commentary*string#organism*string#ligandStructureId*string!..."
or
[{'ecNumber':string, 'inhibitor':string, 'commentary':string, 'organism':string, 'ligandStructureId':string, 'literature':[int1,int2,int3,...]},{'ecNumber':string, 'inhibitor':string, 'commentary':string, 'organism':string, 'ligandStructureId':string, 'literature':[int1,int2,int3,...]},...]"

KCat KM Value
45. getEcNumbersFromKcatKmValue()
Input
top

Input Data Type String (for all languages):
email,password
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password)
resultString = client.service.getEcNumbersFromKcatKmValue(*parameters)
Output
String containing the different EC Numbers (separated by !) or an array of the different EC Numbers (Python 3) linked to KCat KM Value, e.g. "EC Number1!EC Number2!EC Number3!EC Number4!..." or "['EC Number1','EC Number2','EC Number3','EC Number4']"

46. getOrganismsFromKcatKmValue()
Input
top
Input Data Type String (for all languages):
email,password
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password,"ecNumber*", "kcatKmValue*", "kcatKmValueMaximum*", "substrate*", "commentary*", "organism*", "ligandStructureId*", "literature*")
resultString = client.service.getOrganismsFromKcatKmValue(*parameters)
Output
String containing the different Organisms (separated by !) or an array of Organisms (Python 3) linked to KCat KM Value, e.g. "Organism1!Organism2!Organism3!Organism4!..." or "['Organism1','Organism2','Organism3','Organism4']"

47. getKcatKmValue(String)
Input
top
Input Data Type String (for all languages):

email,password,"ecNumber*1.1.1.1#organism*Mus musculus"
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password,"ecNumber*", "kcatKmValue*", "kcatKmValueMaximum*", "substrate*", "commentary*", "organism*", "ligandStructureId*", "literature*")
resultString = client.service.getKcatKmValue(*parameters)
Either the key field ecNumber (e.g. "ecNumber*1.1.1.2"), the key field organism (e.g. "organism*Homo sapiens") or both (e.g. "ecNumber*1.1.1.2#organism*Homo sapiens") have to be specified. If none of these key fields is used as an input parameter, an empty String will be returned by the SOAP query.
In addition, the following optional parameters can be specified: "kcatKmValue", "kcatKmValueMaximum", "substrate", "commentary", "ligandStructureId", "literature
Output
String containing the KCat KM Value entries (entries separated by !, fields separated by #) or an array of KCat KM Value objects (Python 3), e.g. "ecNumber*string#kcatKmValue*string#kcatKmValueMaximum*string#substrate*string#commentary*string#organism*string#ligandStructureId*string!
 ecNumber*string#kcatKmValue*string#kcatKmValueMaximum*string#substrate*string#commentary*string#organism*string#ligandStructureId*string!..."
or
[{'ecNumber':string, 'kcatKmValue':string, 'kcatKmValueMaximum':string, 'substrate':string, 'commentary':string, 'organism':string, 'ligandStructureId':string, 'literature':[int1,int2,int3,...]},{'ecNumber':string, 'kcatKmValue':string, 'kcatKmValueMaximum':string, 'substrate':string, 'commentary':string, 'organism':string, 'ligandStructureId':string, 'literature':[int1,int2,int3,...]},...]"

KI Value
48. getEcNumbersFromKiValue()
Input
top

Input Data Type String (for all languages):
email,password
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password)
resultString = client.service.getEcNumbersFromKiValue(*parameters)
Output
String containing the different EC Numbers (separated by !) or an array of the different EC Numbers (Python 3) linked to KI Value, e.g. "EC Number1!EC Number2!EC Number3!EC Number4!..." or "['EC Number1','EC Number2','EC Number3','EC Number4']"

49. getOrganismsFromKiValue()
Input
top
Input Data Type String (for all languages):
email,password
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password,"ecNumber*", "kiValue*", "kiValueMaximum*", "inhibitor*", "commentary*", "organism*", "ligandStructureId*", "literature*")
resultString = client.service.getOrganismsFromKiValue(*parameters)
Output
String containing the different Organisms (separated by !) or an array of Organisms (Python 3) linked to KI Value, e.g. "Organism1!Organism2!Organism3!Organism4!..." or "['Organism1','Organism2','Organism3','Organism4']"

50. getKiValue(String)
Input
top
Input Data Type String (for all languages):

email,password,"ecNumber*1.1.1.1#organism*Mus musculus"
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password,"ecNumber*", "kiValue*", "kiValueMaximum*", "inhibitor*", "commentary*", "organism*", "ligandStructureId*", "literature*")
resultString = client.service.getKiValue(*parameters)
Either the key field ecNumber (e.g. "ecNumber*1.1.1.2"), the key field organism (e.g. "organism*Homo sapiens") or both (e.g. "ecNumber*1.1.1.2#organism*Homo sapiens") have to be specified. If none of these key fields is used as an input parameter, an empty String will be returned by the SOAP query.
In addition, the following optional parameters can be specified: "kiValue", "kiValueMaximum", "inhibitor", "commentary", "ligandStructureId", "literature
Output
String containing the KI Value entries (entries separated by !, fields separated by #) or an array of KI Value objects (Python 3), e.g. "ecNumber*string#kiValue*string#kiValueMaximum*string#inhibitor*string#commentary*string#organism*string#ligandStructureId*string!
 ecNumber*string#kiValue*string#kiValueMaximum*string#inhibitor*string#commentary*string#organism*string#ligandStructureId*string!..."
or
[{'ecNumber':string, 'kiValue':string, 'kiValueMaximum':string, 'inhibitor':string, 'commentary':string, 'organism':string, 'ligandStructureId':string, 'literature':[int1,int2,int3,...]},{'ecNumber':string, 'kiValue':string, 'kiValueMaximum':string, 'inhibitor':string, 'commentary':string, 'organism':string, 'ligandStructureId':string, 'literature':[int1,int2,int3,...]},...]"

KM Value
51. getEcNumbersFromKmValue()
Input
top

Input Data Type String (for all languages):
email,password
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password)
resultString = client.service.getEcNumbersFromKmValue(*parameters)
Output
String containing the different EC Numbers (separated by !) or an array of the different EC Numbers (Python 3) linked to KM Value, e.g. "EC Number1!EC Number2!EC Number3!EC Number4!..." or "['EC Number1','EC Number2','EC Number3','EC Number4']"

52. getOrganismsFromKmValue()
Input
top
Input Data Type String (for all languages):
email,password
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password,"ecNumber*", "kmValue*", "kmValueMaximum*", "substrate*", "commentary*", "organism*", "ligandStructureId*", "literature*")
resultString = client.service.getOrganismsFromKmValue(*parameters)
Output
String containing the different Organisms (separated by !) or an array of Organisms (Python 3) linked to KM Value, e.g. "Organism1!Organism2!Organism3!Organism4!..." or "['Organism1','Organism2','Organism3','Organism4']"

53. getKmValue(String)
Input
top
Input Data Type String (for all languages):

email,password,"ecNumber*1.1.1.1#organism*Mus musculus"
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password,"ecNumber*", "kmValue*", "kmValueMaximum*", "substrate*", "commentary*", "organism*", "ligandStructureId*", "literature*")
resultString = client.service.getKmValue(*parameters)
Either the key field ecNumber (e.g. "ecNumber*1.1.1.2"), the key field organism (e.g. "organism*Homo sapiens") or both (e.g. "ecNumber*1.1.1.2#organism*Homo sapiens") have to be specified. If none of these key fields is used as an input parameter, an empty String will be returned by the SOAP query.
In addition, the following optional parameters can be specified: "kmValue", "kmValueMaximum", "substrate", "commentary", "ligandStructureId", "literature
Output
String containing the KM Value entries (entries separated by !, fields separated by #) or an array of KM Value objects (Python 3), e.g. "ecNumber*string#kmValue*string#kmValueMaximum*string#substrate*string#commentary*string#organism*string#ligandStructureId*string!
 ecNumber*string#kmValue*string#kmValueMaximum*string#substrate*string#commentary*string#organism*string#ligandStructureId*string!..."
or
[{'ecNumber':string, 'kmValue':string, 'kmValueMaximum':string, 'substrate':string, 'commentary':string, 'organism':string, 'ligandStructureId':string, 'literature':[int1,int2,int3,...]},{'ecNumber':string, 'kmValue':string, 'kmValueMaximum':string, 'substrate':string, 'commentary':string, 'organism':string, 'ligandStructureId':string, 'literature':[int1,int2,int3,...]},...]"

Ligands
54. getEcNumbersFromLigands()
Input
top

Input Data Type String (for all languages):
email,password
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password)
resultString = client.service.getEcNumbersFromLigands(*parameters)
Output
String containing the different EC Numbers (separated by !) or an array of the different EC Numbers (Python 3) linked to Ligands, e.g. "EC Number1!EC Number2!EC Number3!EC Number4!..." or "['EC Number1','EC Number2','EC Number3','EC Number4']"

55. getOrganismsFromLigands()
Input
top
Input Data Type String (for all languages):
email,password
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password,"ecNumber*", "role*", "ligand*", "organism*", "ligandStructureId*")
resultString = client.service.getOrganismsFromLigands(*parameters)
Output
String containing the different Organisms (separated by !) or an array of Organisms (Python 3) linked to Ligands, e.g. "Organism1!Organism2!Organism3!Organism4!..." or "['Organism1','Organism2','Organism3','Organism4']"

56. getLigands(String)
Input
top
Input Data Type String (for all languages):

email,password,"ecNumber*1.1.1.1#organism*Mus musculus"
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password,"ecNumber*", "role*", "ligand*", "organism*", "ligandStructureId*")
resultString = client.service.getLigands(*parameters)
Either the key field ecNumber (e.g. "ecNumber*1.1.1.2"), the key field organism (e.g. "organism*Homo sapiens") or both (e.g. "ecNumber*1.1.1.2#organism*Homo sapiens") have to be specified. If none of these key fields is used as an input parameter, an empty String will be returned by the SOAP query.
In addition, the following optional parameters can be specified: "role", "ligand", "ligandStructureId
Output
String containing the Ligands entries (entries separated by !, fields separated by #) or an array of Ligands objects (Python 3), e.g. "ecNumber*string#role*string#ligand*string#organism*string#ligandStructureId*string!
 ecNumber*string#role*string#ligand*string#organism*string#ligandStructureId*string!..."
or
[{'ecNumber':string, 'role':string, 'ligand':string, 'organism':string, 'ligandStructureId':string},{'ecNumber':string, 'role':string, 'ligand':string, 'organism':string, 'ligandStructureId':string},...]"

Localization
57. getEcNumbersFromLocalization()
Input
top

Input Data Type String (for all languages):
email,password
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password)
resultString = client.service.getEcNumbersFromLocalization(*parameters)
Output
String containing the different EC Numbers (separated by !) or an array of the different EC Numbers (Python 3) linked to Localization, e.g. "EC Number1!EC Number2!EC Number3!EC Number4!..." or "['EC Number1','EC Number2','EC Number3','EC Number4']"

58. getOrganismsFromLocalization()
Input
top
Input Data Type String (for all languages):
email,password
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password,"ecNumber*", "localization*", "commentary*", "organism*", "idGo*", "literature*", "textmining*")
resultString = client.service.getOrganismsFromLocalization(*parameters)
Output
String containing the different Organisms (separated by !) or an array of Organisms (Python 3) linked to Localization, e.g. "Organism1!Organism2!Organism3!Organism4!..." or "['Organism1','Organism2','Organism3','Organism4']"

59. getLocalization(String)
Input
top
Input Data Type String (for all languages):

email,password,"ecNumber*1.1.1.1#organism*Mus musculus"
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password,"ecNumber*", "localization*", "commentary*", "organism*", "idGo*", "literature*", "textmining*")
resultString = client.service.getLocalization(*parameters)
Either the key field ecNumber (e.g. "ecNumber*1.1.1.2"), the key field organism (e.g. "organism*Homo sapiens") or both (e.g. "ecNumber*1.1.1.2#organism*Homo sapiens") have to be specified. If none of these key fields is used as an input parameter, an empty String will be returned by the SOAP query.
In addition, the following optional parameters can be specified: "localization", "commentary", "idGo", "literature", "textmining
Setting the field/parameter textmining to zero ("textmining*0"), restricts the search to manually annotated entries only.
Output
String containing the Localization entries (entries separated by !, fields separated by #) or an array of Localization objects (Python 3), e.g. "ecNumber*string#localization*string#commentary*string#organism*string#idGo*string#textmining*string!
 ecNumber*string#localization*string#commentary*string#organism*string#idGo*string#textmining*string!..."
or
[{'ecNumber':string, 'localization':string, 'commentary':string, 'organism':string, 'idGo':string, 'literature':[int1,int2,int3,...], 'textmining':string},{'ecNumber':string, 'localization':string, 'commentary':string, 'organism':string, 'idGo':string, 'literature':[int1,int2,int3,...], 'textmining':string},...]"

Metals Ions
60. getEcNumbersFromMetalsIons()
Input
top

Input Data Type String (for all languages):
email,password
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password)
resultString = client.service.getEcNumbersFromMetalsIons(*parameters)
Output
String containing the different EC Numbers (separated by !) or an array of the different EC Numbers (Python 3) linked to Metals Ions, e.g. "EC Number1!EC Number2!EC Number3!EC Number4!..." or "['EC Number1','EC Number2','EC Number3','EC Number4']"

61. getOrganismsFromMetalsIons()
Input
top
Input Data Type String (for all languages):
email,password
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password,"ecNumber*", "metalsIons*", "commentary*", "organism*", "ligandStructureId*", "literature*")
resultString = client.service.getOrganismsFromMetalsIons(*parameters)
Output
String containing the different Organisms (separated by !) or an array of Organisms (Python 3) linked to Metals Ions, e.g. "Organism1!Organism2!Organism3!Organism4!..." or "['Organism1','Organism2','Organism3','Organism4']"

62. getMetalsIons(String)
Input
top
Input Data Type String (for all languages):

email,password,"ecNumber*1.1.1.1#organism*Mus musculus"
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password,"ecNumber*", "metalsIons*", "commentary*", "organism*", "ligandStructureId*", "literature*")
resultString = client.service.getMetalsIons(*parameters)
Either the key field ecNumber (e.g. "ecNumber*1.1.1.2"), the key field organism (e.g. "organism*Homo sapiens") or both (e.g. "ecNumber*1.1.1.2#organism*Homo sapiens") have to be specified. If none of these key fields is used as an input parameter, an empty String will be returned by the SOAP query.
In addition, the following optional parameters can be specified: "metalsIons", "commentary", "ligandStructureId", "literature
Output
String containing the Metals Ions entries (entries separated by !, fields separated by #) or an array of Metals Ions objects (Python 3), e.g. "ecNumber*string#metalsIons*string#commentary*string#organism*string#ligandStructureId*string!
 ecNumber*string#metalsIons*string#commentary*string#organism*string#ligandStructureId*string!..."
or
[{'ecNumber':string, 'metalsIons':string, 'commentary':string, 'organism':string, 'ligandStructureId':string, 'literature':[int1,int2,int3,...]},{'ecNumber':string, 'metalsIons':string, 'commentary':string, 'organism':string, 'ligandStructureId':string, 'literature':[int1,int2,int3,...]},...]"

Molecular Weight
63. getEcNumbersFromMolecularWeight()
Input
top

Input Data Type String (for all languages):
email,password
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password)
resultString = client.service.getEcNumbersFromMolecularWeight(*parameters)
Output
String containing the different EC Numbers (separated by !) or an array of the different EC Numbers (Python 3) linked to Molecular Weight, e.g. "EC Number1!EC Number2!EC Number3!EC Number4!..." or "['EC Number1','EC Number2','EC Number3','EC Number4']"

64. getOrganismsFromMolecularWeight()
Input
top
Input Data Type String (for all languages):
email,password
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password,"ecNumber*", "molecularWeight*", "molecularWeightMaximum*", "commentary*", "organism*", "literature*")
resultString = client.service.getOrganismsFromMolecularWeight(*parameters)
Output
String containing the different Organisms (separated by !) or an array of Organisms (Python 3) linked to Molecular Weight, e.g. "Organism1!Organism2!Organism3!Organism4!..." or "['Organism1','Organism2','Organism3','Organism4']"

65. getMolecularWeight(String)
Input
top
Input Data Type String (for all languages):

email,password,"ecNumber*1.1.1.1#organism*Mus musculus"
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password,"ecNumber*", "molecularWeight*", "molecularWeightMaximum*", "commentary*", "organism*", "literature*")
resultString = client.service.getMolecularWeight(*parameters)
Either the key field ecNumber (e.g. "ecNumber*1.1.1.2"), the key field organism (e.g. "organism*Homo sapiens") or both (e.g. "ecNumber*1.1.1.2#organism*Homo sapiens") have to be specified. If none of these key fields is used as an input parameter, an empty String will be returned by the SOAP query.
In addition, the following optional parameters can be specified: "molecularWeight", "molecularWeightMaximum", "commentary", "literature
Output
String containing the Molecular Weight entries (entries separated by !, fields separated by #) or an array of Molecular Weight objects (Python 3), e.g. "ecNumber*string#molecularWeight*string#molecularWeightMaximum*string#commentary*string#organism*string!
 ecNumber*string#molecularWeight*string#molecularWeightMaximum*string#commentary*string#organism*string!..."
or
[{'ecNumber':string, 'molecularWeight':string, 'molecularWeightMaximum':string, 'commentary':string, 'organism':string, 'literature':[int1,int2,int3,...]},{'ecNumber':string, 'molecularWeight':string, 'molecularWeightMaximum':string, 'commentary':string, 'organism':string, 'literature':[int1,int2,int3,...]},...]"

Natural Product
66. getEcNumbersFromNaturalProduct()
Input
top

Input Data Type String (for all languages):
email,password
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password)
resultString = client.service.getEcNumbersFromNaturalProduct(*parameters)
Output
String containing the different EC Numbers (separated by !) or an array of the different EC Numbers (Python 3) linked to Natural Product, e.g. "EC Number1!EC Number2!EC Number3!EC Number4!..." or "['EC Number1','EC Number2','EC Number3','EC Number4']"

67. getOrganismsFromNaturalProduct()
Input
top
Input Data Type String (for all languages):
email,password
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password,"ecNumber*", "naturalProduct*", "naturalReactionPartners*", "organism*", "ligandStructureId*")
resultString = client.service.getOrganismsFromNaturalProduct(*parameters)
Output
String containing the different Organisms (separated by !) or an array of Organisms (Python 3) linked to Natural Product, e.g. "Organism1!Organism2!Organism3!Organism4!..." or "['Organism1','Organism2','Organism3','Organism4']"

68. getNaturalProduct(String)
Input
top
Input Data Type String (for all languages):

email,password,"ecNumber*1.1.1.1#organism*Mus musculus"
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password,"ecNumber*", "naturalProduct*", "naturalReactionPartners*", "organism*", "ligandStructureId*")
resultString = client.service.getNaturalProduct(*parameters)
Either the key field ecNumber (e.g. "ecNumber*1.1.1.2"), the key field organism (e.g. "organism*Homo sapiens") or both (e.g. "ecNumber*1.1.1.2#organism*Homo sapiens") have to be specified. If none of these key fields is used as an input parameter, an empty String will be returned by the SOAP query.
In addition, the following optional parameters can be specified: "naturalProduct", "naturalReactionPartners", "ligandStructureId
Output
String containing the Natural Product entries (entries separated by !, fields separated by #) or an array of Natural Product objects (Python 3), e.g. "ecNumber*string#naturalProduct*string#naturalReactionPartners*string#organism*string#ligandStructureId*string!
 ecNumber*string#naturalProduct*string#naturalReactionPartners*string#organism*string#ligandStructureId*string!..."
or
[{'ecNumber':string, 'naturalProduct':string, 'naturalReactionPartners':string, 'organism':string, 'ligandStructureId':string},{'ecNumber':string, 'naturalProduct':string, 'naturalReactionPartners':string, 'organism':string, 'ligandStructureId':string},...]"

Natural Substrate
69. getEcNumbersFromNaturalSubstrate()
Input
top

Input Data Type String (for all languages):
email,password
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password)
resultString = client.service.getEcNumbersFromNaturalSubstrate(*parameters)
Output
String containing the different EC Numbers (separated by !) or an array of the different EC Numbers (Python 3) linked to Natural Substrate, e.g. "EC Number1!EC Number2!EC Number3!EC Number4!..." or "['EC Number1','EC Number2','EC Number3','EC Number4']"

70. getOrganismsFromNaturalSubstrate()
Input
top
Input Data Type String (for all languages):
email,password
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password,"ecNumber*", "naturalSubstrate*", "naturalReactionPartners*", "organism*", "ligandStructureId*")
resultString = client.service.getOrganismsFromNaturalSubstrate(*parameters)
Output
String containing the different Organisms (separated by !) or an array of Organisms (Python 3) linked to Natural Substrate, e.g. "Organism1!Organism2!Organism3!Organism4!..." or "['Organism1','Organism2','Organism3','Organism4']"

71. getNaturalSubstrate(String)
Input
top
Input Data Type String (for all languages):

email,password,"ecNumber*1.1.1.1#organism*Mus musculus"
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password,"ecNumber*", "naturalSubstrate*", "naturalReactionPartners*", "organism*", "ligandStructureId*")
resultString = client.service.getNaturalSubstrate(*parameters)
Either the key field ecNumber (e.g. "ecNumber*1.1.1.2"), the key field organism (e.g. "organism*Homo sapiens") or both (e.g. "ecNumber*1.1.1.2#organism*Homo sapiens") have to be specified. If none of these key fields is used as an input parameter, an empty String will be returned by the SOAP query.
In addition, the following optional parameters can be specified: "naturalSubstrate", "naturalReactionPartners", "ligandStructureId
Output
String containing the Natural Substrate entries (entries separated by !, fields separated by #) or an array of Natural Substrate objects (Python 3), e.g. "ecNumber*string#naturalSubstrate*string#naturalReactionPartners*string#organism*string#ligandStructureId*string!
 ecNumber*string#naturalSubstrate*string#naturalReactionPartners*string#organism*string#ligandStructureId*string!..."
or
[{'ecNumber':string, 'naturalSubstrate':string, 'naturalReactionPartners':string, 'organism':string, 'ligandStructureId':string},{'ecNumber':string, 'naturalSubstrate':string, 'naturalReactionPartners':string, 'organism':string, 'ligandStructureId':string},...]"

Natural Substrates Products
72. getEcNumbersFromNaturalSubstratesProducts()
Input
top

Input Data Type String (for all languages):
email,password
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password)
resultString = client.service.getEcNumbersFromNaturalSubstratesProducts(*parameters)
Output
String containing the different EC Numbers (separated by !) or an array of the different EC Numbers (Python 3) linked to Natural Substrates Products, e.g. "EC Number1!EC Number2!EC Number3!EC Number4!..." or "['EC Number1','EC Number2','EC Number3','EC Number4']"

73. getNaturalSubstratesProducts(String)
Input
top
Input Data Type String (for all languages):

email,password,"ecNumber*1.1.1.1"
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password,"ecNumber*", "naturalSubstrates*", "organismNaturalSubstrates*", "commentaryNaturalSubstrates*", "naturalProducts*", "commentaryNaturalProducts*", "organismNaturalProducts*", "reversibility*")
resultString = client.service.getNaturalSubstratesProducts(*parameters)
Either the key field ecNumber (e.g. "ecNumber*1.1.1.2"), the key field organism (e.g. "organism*Homo sapiens") or both (e.g. "ecNumber*1.1.1.2#organism*Homo sapiens") have to be specified. If none of these key fields is used as an input parameter, an empty String will be returned by the SOAP query.
In addition, the following optional parameters can be specified: "naturalSubstrates", "organismNaturalSubstrates", "commentaryNaturalSubstrates", "naturalProducts", "commentaryNaturalProducts", "organismNaturalProducts", "reversibility
Output
String containing the Natural Substrates Products entries (entries separated by !, fields separated by #) or an array of Natural Substrates Products objects (Python 3), e.g. "ecNumber*string#naturalSubstrates*string#organismNaturalSubstrates*string#commentaryNaturalSubstrates*string#naturalProducts*string#commentaryNaturalProducts*string#organismNaturalProducts*string#reversibility*string!
 ecNumber*string#naturalSubstrates*string#organismNaturalSubstrates*string#commentaryNaturalSubstrates*string#naturalProducts*string#commentaryNaturalProducts*string#organismNaturalProducts*string#reversibility*string!..."
or
[{'ecNumber':string, 'naturalSubstrates':string, 'organismNaturalSubstrates':string, 'commentaryNaturalSubstrates':string, 'naturalProducts':string, 'commentaryNaturalProducts':string, 'organismNaturalProducts':string, 'reversibility':string},{'ecNumber':string, 'naturalSubstrates':string, 'organismNaturalSubstrates':string, 'commentaryNaturalSubstrates':string, 'naturalProducts':string, 'commentaryNaturalProducts':string, 'organismNaturalProducts':string, 'reversibility':string},...]"

Organic Solvent Stability
74. getEcNumbersFromOrganicSolventStability()
Input
top

Input Data Type String (for all languages):
email,password
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password)
resultString = client.service.getEcNumbersFromOrganicSolventStability(*parameters)
Output
String containing the different EC Numbers (separated by !) or an array of the different EC Numbers (Python 3) linked to Organic Solvent Stability, e.g. "EC Number1!EC Number2!EC Number3!EC Number4!..." or "['EC Number1','EC Number2','EC Number3','EC Number4']"

75. getOrganismsFromOrganicSolventStability()
Input
top
Input Data Type String (for all languages):
email,password
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password,"ecNumber*", "organicSolvent*", "commentary*", "organism*", "literature*")
resultString = client.service.getOrganismsFromOrganicSolventStability(*parameters)
Output
String containing the different Organisms (separated by !) or an array of Organisms (Python 3) linked to Organic Solvent Stability, e.g. "Organism1!Organism2!Organism3!Organism4!..." or "['Organism1','Organism2','Organism3','Organism4']"

76. getOrganicSolventStability(String)
Input
top
Input Data Type String (for all languages):

email,password,"ecNumber*1.1.1.1#organism*Mus musculus"
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password,"ecNumber*", "organicSolvent*", "commentary*", "organism*", "literature*")
resultString = client.service.getOrganicSolventStability(*parameters)
Either the key field ecNumber (e.g. "ecNumber*1.1.1.2"), the key field organism (e.g. "organism*Homo sapiens") or both (e.g. "ecNumber*1.1.1.2#organism*Homo sapiens") have to be specified. If none of these key fields is used as an input parameter, an empty String will be returned by the SOAP query.
In addition, the following optional parameters can be specified: "organicSolvent", "commentary", "literature
Output
String containing the Organic Solvent Stability entries (entries separated by !, fields separated by #) or an array of Organic Solvent Stability objects (Python 3), e.g. "ecNumber*string#organicSolvent*string#commentary*string#organism*string!
 ecNumber*string#organicSolvent*string#commentary*string#organism*string!..."
or
[{'ecNumber':string, 'organicSolvent':string, 'commentary':string, 'organism':string, 'literature':[int1,int2,int3,...]},{'ecNumber':string, 'organicSolvent':string, 'commentary':string, 'organism':string, 'literature':[int1,int2,int3,...]},...]"

Organism
77. getEcNumbersFromOrganism()
Input
top

Input Data Type String (for all languages):
email,password
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password)
resultString = client.service.getEcNumbersFromOrganism(*parameters)
Output
String containing the different EC Numbers (separated by !) or an array of the different EC Numbers (Python 3) linked to Organism, e.g. "EC Number1!EC Number2!EC Number3!EC Number4!..." or "['EC Number1','EC Number2','EC Number3','EC Number4']"

78. getOrganismsFromOrganism()
Input
top
Input Data Type String (for all languages):
email,password
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password,"organism*", "commentary*", "literature*", "sequenceCode*", "textmining*", "ecNumber*")
resultString = client.service.getOrganismsFromOrganism(*parameters)
Output
String containing the different Organisms (separated by !) or an array of Organisms (Python 3) linked to Organism, e.g. "Organism1!Organism2!Organism3!Organism4!..." or "['Organism1','Organism2','Organism3','Organism4']"

79. getOrganism(String)
Input
top
Input Data Type String (for all languages):

email,password,"ecNumber*1.1.1.1#organism*Mus musculus"
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password,"organism*", "commentary*", "literature*", "sequenceCode*", "textmining*", "ecNumber*")
resultString = client.service.getOrganism(*parameters)
Either the key field ecNumber (e.g. "ecNumber*1.1.1.2"), the key field organism (e.g. "organism*Homo sapiens") or both (e.g. "ecNumber*1.1.1.2#organism*Homo sapiens") have to be specified. If none of these key fields is used as an input parameter, an empty String will be returned by the SOAP query.
In addition, the following optional parameters can be specified: "commentary", "literature", "sequenceCode", "textmining
Setting the field/parameter textmining to zero ("textmining*0"), restricts the search to manually annotated entries only.
Output
String containing the Organism entries (entries separated by !, fields separated by #) or an array of Organism objects (Python 3), e.g. "organism*string#commentary*string#sequenceCode*string#textmining*string#ecNumber*string!
 organism*string#commentary*string#sequenceCode*string#textmining*string#ecNumber*string!..."
or
[{'organism':string, 'commentary':string, 'literature':[int1,int2,int3,...], 'sequenceCode':string, 'textmining':string, 'ecNumber':string},{'organism':string, 'commentary':string, 'literature':[int1,int2,int3,...], 'sequenceCode':string, 'textmining':string, 'ecNumber':string},...]"

Oxidation Stability
80. getEcNumbersFromOxidationStability()
Input
top

Input Data Type String (for all languages):
email,password
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password)
resultString = client.service.getEcNumbersFromOxidationStability(*parameters)
Output
String containing the different EC Numbers (separated by !) or an array of the different EC Numbers (Python 3) linked to Oxidation Stability, e.g. "EC Number1!EC Number2!EC Number3!EC Number4!..." or "['EC Number1','EC Number2','EC Number3','EC Number4']"

81. getOrganismsFromOxidationStability()
Input
top
Input Data Type String (for all languages):
email,password
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password,"ecNumber*", "oxidationStability*", "organism*", "literature*")
resultString = client.service.getOrganismsFromOxidationStability(*parameters)
Output
String containing the different Organisms (separated by !) or an array of Organisms (Python 3) linked to Oxidation Stability, e.g. "Organism1!Organism2!Organism3!Organism4!..." or "['Organism1','Organism2','Organism3','Organism4']"

82. getOxidationStability(String)
Input
top
Input Data Type String (for all languages):

email,password,"ecNumber*1.1.1.1#organism*Mus musculus"
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password,"ecNumber*", "oxidationStability*", "organism*", "literature*")
resultString = client.service.getOxidationStability(*parameters)
Either the key field ecNumber (e.g. "ecNumber*1.1.1.2"), the key field organism (e.g. "organism*Homo sapiens") or both (e.g. "ecNumber*1.1.1.2#organism*Homo sapiens") have to be specified. If none of these key fields is used as an input parameter, an empty String will be returned by the SOAP query.
In addition, the following optional parameters can be specified: "oxidationStability", "literature
Output
String containing the Oxidation Stability entries (entries separated by !, fields separated by #) or an array of Oxidation Stability objects (Python 3), e.g. "ecNumber*string#oxidationStability*string#organism*string!
 ecNumber*string#oxidationStability*string#organism*string!..."
or
[{'ecNumber':string, 'oxidationStability':string, 'organism':string, 'literature':[int1,int2,int3,...]},{'ecNumber':string, 'oxidationStability':string, 'organism':string, 'literature':[int1,int2,int3,...]},...]"

Pathway
83. getEcNumbersFromPathway()
Input
top

Input Data Type String (for all languages):
email,password
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password)
resultString = client.service.getEcNumbersFromPathway(*parameters)
Output
String containing the different EC Numbers (separated by !) or an array of the different EC Numbers (Python 3) linked to Pathway, e.g. "EC Number1!EC Number2!EC Number3!EC Number4!..." or "['EC Number1','EC Number2','EC Number3','EC Number4']"

84. getPathway(String)
Input
top
Input Data Type String (for all languages):

email,password,"ecNumber*1.1.1.1"
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password,"ecNumber*", "pathway*", "link*", "source_database*")
resultString = client.service.getPathway(*parameters)
Either the key field ecNumber (e.g. "ecNumber*1.1.1.2"), the key field organism (e.g. "organism*Homo sapiens") or both (e.g. "ecNumber*1.1.1.2#organism*Homo sapiens") have to be specified. If none of these key fields is used as an input parameter, an empty String will be returned by the SOAP query.
In addition, the following optional parameters can be specified: "pathway", "link", "source_database
Output
String containing the Pathway entries (entries separated by !, fields separated by #) or an array of Pathway objects (Python 3), e.g. "ecNumber*string#pathway*string#link*string#source_database*string!
 ecNumber*string#pathway*string#link*string#source_database*string!..."
or
[{'ecNumber':string, 'pathway':string, 'link':string, 'source_database':string},{'ecNumber':string, 'pathway':string, 'link':string, 'source_database':string},...]"

PDB
85. getEcNumbersFromPdb()
Input
top

Input Data Type String (for all languages):
email,password
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password)
resultString = client.service.getEcNumbersFromPdb(*parameters)
Output
String containing the different EC Numbers (separated by !) or an array of the different EC Numbers (Python 3) linked to PDB, e.g. "EC Number1!EC Number2!EC Number3!EC Number4!..." or "['EC Number1','EC Number2','EC Number3','EC Number4']"

86. getOrganismsFromPdb()
Input
top
Input Data Type String (for all languages):
email,password
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password,"ecNumber*", "pdb*", "organism*")
resultString = client.service.getOrganismsFromPdb(*parameters)
Output
String containing the different Organisms (separated by !) or an array of Organisms (Python 3) linked to PDB, e.g. "Organism1!Organism2!Organism3!Organism4!..." or "['Organism1','Organism2','Organism3','Organism4']"

87. getPdb(String)
Input
top
Input Data Type String (for all languages):

email,password,"ecNumber*1.1.1.1#organism*Mus musculus"
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password,"ecNumber*", "pdb*", "organism*")
resultString = client.service.getPdb(*parameters)
Either the key field ecNumber (e.g. "ecNumber*1.1.1.2"), the key field organism (e.g. "organism*Homo sapiens") or both (e.g. "ecNumber*1.1.1.2#organism*Homo sapiens") have to be specified. If none of these key fields is used as an input parameter, an empty String will be returned by the SOAP query.
In addition, the following optional parameters can be specified: "pdb
Output
String containing the PDB entries (entries separated by !, fields separated by #) or an array of PDB objects (Python 3), e.g. "ecNumber*string#pdb*string#organism*string!
 ecNumber*string#pdb*string#organism*string!..."
or
[{'ecNumber':string, 'pdb':string, 'organism':string},{'ecNumber':string, 'pdb':string, 'organism':string},...]"

pH Optimum
88. getEcNumbersFromPhOptimum()
Input
top

Input Data Type String (for all languages):
email,password
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password)
resultString = client.service.getEcNumbersFromPhOptimum(*parameters)
Output
String containing the different EC Numbers (separated by !) or an array of the different EC Numbers (Python 3) linked to pH Optimum, e.g. "EC Number1!EC Number2!EC Number3!EC Number4!..." or "['EC Number1','EC Number2','EC Number3','EC Number4']"

89. getOrganismsFromPhOptimum()
Input
top
Input Data Type String (for all languages):
email,password
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password,"ecNumber*", "phOptimum*", "phOptimumMaximum*", "commentary*", "organism*", "literature*")
resultString = client.service.getOrganismsFromPhOptimum(*parameters)
Output
String containing the different Organisms (separated by !) or an array of Organisms (Python 3) linked to pH Optimum, e.g. "Organism1!Organism2!Organism3!Organism4!..." or "['Organism1','Organism2','Organism3','Organism4']"

90. getPhOptimum(String)
Input
top
Input Data Type String (for all languages):

email,password,"ecNumber*1.1.1.1#organism*Mus musculus"
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password,"ecNumber*", "phOptimum*", "phOptimumMaximum*", "commentary*", "organism*", "literature*")
resultString = client.service.getPhOptimum(*parameters)
Either the key field ecNumber (e.g. "ecNumber*1.1.1.2"), the key field organism (e.g. "organism*Homo sapiens") or both (e.g. "ecNumber*1.1.1.2#organism*Homo sapiens") have to be specified. If none of these key fields is used as an input parameter, an empty String will be returned by the SOAP query.
In addition, the following optional parameters can be specified: "phOptimum", "phOptimumMaximum", "commentary", "literature
Output
String containing the pH Optimum entries (entries separated by !, fields separated by #) or an array of pH Optimum objects (Python 3), e.g. "ecNumber*string#phOptimum*string#phOptimumMaximum*string#commentary*string#organism*string!
 ecNumber*string#phOptimum*string#phOptimumMaximum*string#commentary*string#organism*string!..."
or
[{'ecNumber':string, 'phOptimum':string, 'phOptimumMaximum':string, 'commentary':string, 'organism':string, 'literature':[int1,int2,int3,...]},{'ecNumber':string, 'phOptimum':string, 'phOptimumMaximum':string, 'commentary':string, 'organism':string, 'literature':[int1,int2,int3,...]},...]"

pH Range
91. getEcNumbersFromPhRange()
Input
top

Input Data Type String (for all languages):
email,password
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password)
resultString = client.service.getEcNumbersFromPhRange(*parameters)
Output
String containing the different EC Numbers (separated by !) or an array of the different EC Numbers (Python 3) linked to pH Range, e.g. "EC Number1!EC Number2!EC Number3!EC Number4!..." or "['EC Number1','EC Number2','EC Number3','EC Number4']"

92. getOrganismsFromPhRange()
Input
top
Input Data Type String (for all languages):
email,password
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password,"ecNumber*", "phRange*", "phRangeMaximum*", "commentary*", "organism*", "literature*")
resultString = client.service.getOrganismsFromPhRange(*parameters)
Output
String containing the different Organisms (separated by !) or an array of Organisms (Python 3) linked to pH Range, e.g. "Organism1!Organism2!Organism3!Organism4!..." or "['Organism1','Organism2','Organism3','Organism4']"

93. getPhRange(String)
Input
top
Input Data Type String (for all languages):

email,password,"ecNumber*1.1.1.1#organism*Mus musculus"
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password,"ecNumber*", "phRange*", "phRangeMaximum*", "commentary*", "organism*", "literature*")
resultString = client.service.getPhRange(*parameters)
Either the key field ecNumber (e.g. "ecNumber*1.1.1.2"), the key field organism (e.g. "organism*Homo sapiens") or both (e.g. "ecNumber*1.1.1.2#organism*Homo sapiens") have to be specified. If none of these key fields is used as an input parameter, an empty String will be returned by the SOAP query.
In addition, the following optional parameters can be specified: "phRange", "phRangeMaximum", "commentary", "literature
Output
String containing the pH Range entries (entries separated by !, fields separated by #) or an array of pH Range objects (Python 3), e.g. "ecNumber*string#phRange*string#phRangeMaximum*string#commentary*string#organism*string!
 ecNumber*string#phRange*string#phRangeMaximum*string#commentary*string#organism*string!..."
or
[{'ecNumber':string, 'phRange':string, 'phRangeMaximum':string, 'commentary':string, 'organism':string, 'literature':[int1,int2,int3,...]},{'ecNumber':string, 'phRange':string, 'phRangeMaximum':string, 'commentary':string, 'organism':string, 'literature':[int1,int2,int3,...]},...]"

pH Stability
94. getEcNumbersFromPhStability()
Input
top

Input Data Type String (for all languages):
email,password
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password)
resultString = client.service.getEcNumbersFromPhStability(*parameters)
Output
String containing the different EC Numbers (separated by !) or an array of the different EC Numbers (Python 3) linked to pH Stability, e.g. "EC Number1!EC Number2!EC Number3!EC Number4!..." or "['EC Number1','EC Number2','EC Number3','EC Number4']"

95. getOrganismsFromPhStability()
Input
top
Input Data Type String (for all languages):
email,password
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password,"ecNumber*", "phStability*", "phStabilityMaximum*", "commentary*", "organism*", "literature*")
resultString = client.service.getOrganismsFromPhStability(*parameters)
Output
String containing the different Organisms (separated by !) or an array of Organisms (Python 3) linked to pH Stability, e.g. "Organism1!Organism2!Organism3!Organism4!..." or "['Organism1','Organism2','Organism3','Organism4']"

96. getPhStability(String)
Input
top
Input Data Type String (for all languages):

email,password,"ecNumber*1.1.1.1#organism*Mus musculus"
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password,"ecNumber*", "phStability*", "phStabilityMaximum*", "commentary*", "organism*", "literature*")
resultString = client.service.getPhStability(*parameters)
Either the key field ecNumber (e.g. "ecNumber*1.1.1.2"), the key field organism (e.g. "organism*Homo sapiens") or both (e.g. "ecNumber*1.1.1.2#organism*Homo sapiens") have to be specified. If none of these key fields is used as an input parameter, an empty String will be returned by the SOAP query.
In addition, the following optional parameters can be specified: "phStability", "phStabilityMaximum", "commentary", "literature
Output
String containing the pH Stability entries (entries separated by !, fields separated by #) or an array of pH Stability objects (Python 3), e.g. "ecNumber*string#phStability*string#phStabilityMaximum*string#commentary*string#organism*string!
 ecNumber*string#phStability*string#phStabilityMaximum*string#commentary*string#organism*string!..."
or
[{'ecNumber':string, 'phStability':string, 'phStabilityMaximum':string, 'commentary':string, 'organism':string, 'literature':[int1,int2,int3,...]},{'ecNumber':string, 'phStability':string, 'phStabilityMaximum':string, 'commentary':string, 'organism':string, 'literature':[int1,int2,int3,...]},...]"

pI Value
97. getEcNumbersFromPiValue()
Input
top

Input Data Type String (for all languages):
email,password
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password)
resultString = client.service.getEcNumbersFromPiValue(*parameters)
Output
String containing the different EC Numbers (separated by !) or an array of the different EC Numbers (Python 3) linked to pI Value, e.g. "EC Number1!EC Number2!EC Number3!EC Number4!..." or "['EC Number1','EC Number2','EC Number3','EC Number4']"

98. getOrganismsFromPiValue()
Input
top
Input Data Type String (for all languages):
email,password
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password,"ecNumber*", "piValue*", "piValueMaximum*", "commentary*", "literature*", "organism*")
resultString = client.service.getOrganismsFromPiValue(*parameters)
Output
String containing the different Organisms (separated by !) or an array of Organisms (Python 3) linked to pI Value, e.g. "Organism1!Organism2!Organism3!Organism4!..." or "['Organism1','Organism2','Organism3','Organism4']"

99. getPiValue(String)
Input
top
Input Data Type String (for all languages):

email,password,"ecNumber*1.1.1.1#organism*Mus musculus"
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password,"ecNumber*", "piValue*", "piValueMaximum*", "commentary*", "literature*", "organism*")
resultString = client.service.getPiValue(*parameters)
Either the key field ecNumber (e.g. "ecNumber*1.1.1.2"), the key field organism (e.g. "organism*Homo sapiens") or both (e.g. "ecNumber*1.1.1.2#organism*Homo sapiens") have to be specified. If none of these key fields is used as an input parameter, an empty String will be returned by the SOAP query.
In addition, the following optional parameters can be specified: "piValue", "piValueMaximum", "commentary", "literature
Output
String containing the pI Value entries (entries separated by !, fields separated by #) or an array of pI Value objects (Python 3), e.g. "ecNumber*string#piValue*string#piValueMaximum*string#commentary*string#organism*string!
 ecNumber*string#piValue*string#piValueMaximum*string#commentary*string#organism*string!..."
or
[{'ecNumber':string, 'piValue':string, 'piValueMaximum':string, 'commentary':string, 'literature':[int1,int2,int3,...], 'organism':string},{'ecNumber':string, 'piValue':string, 'piValueMaximum':string, 'commentary':string, 'literature':[int1,int2,int3,...], 'organism':string},...]"

Posttranslational Modification
100. getEcNumbersFromPosttranslationalModification()
Input
top

Input Data Type String (for all languages):
email,password
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password)
resultString = client.service.getEcNumbersFromPosttranslationalModification(*parameters)
Output
String containing the different EC Numbers (separated by !) or an array of the different EC Numbers (Python 3) linked to Posttranslational Modification, e.g. "EC Number1!EC Number2!EC Number3!EC Number4!..." or "['EC Number1','EC Number2','EC Number3','EC Number4']"

101. getOrganismsFromPosttranslationalModification()
Input
top
Input Data Type String (for all languages):
email,password
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password,"ecNumber*", "posttranslationalModification*", "commentary*", "organism*", "literature*")
resultString = client.service.getOrganismsFromPosttranslationalModification(*parameters)
Output
String containing the different Organisms (separated by !) or an array of Organisms (Python 3) linked to Posttranslational Modification, e.g. "Organism1!Organism2!Organism3!Organism4!..." or "['Organism1','Organism2','Organism3','Organism4']"

102. getPosttranslationalModification(String)
Input
top
Input Data Type String (for all languages):

email,password,"ecNumber*1.1.1.1#organism*Mus musculus"
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password,"ecNumber*", "posttranslationalModification*", "commentary*", "organism*", "literature*")
resultString = client.service.getPosttranslationalModification(*parameters)
Either the key field ecNumber (e.g. "ecNumber*1.1.1.2"), the key field organism (e.g. "organism*Homo sapiens") or both (e.g. "ecNumber*1.1.1.2#organism*Homo sapiens") have to be specified. If none of these key fields is used as an input parameter, an empty String will be returned by the SOAP query.
In addition, the following optional parameters can be specified: "posttranslationalModification", "commentary", "literature
Output
String containing the Posttranslational Modification entries (entries separated by !, fields separated by #) or an array of Posttranslational Modification objects (Python 3), e.g. "ecNumber*string#posttranslationalModification*string#commentary*string#organism*string!
 ecNumber*string#posttranslationalModification*string#commentary*string#organism*string!..."
or
[{'ecNumber':string, 'posttranslationalModification':string, 'commentary':string, 'organism':string, 'literature':[int1,int2,int3,...]},{'ecNumber':string, 'posttranslationalModification':string, 'commentary':string, 'organism':string, 'literature':[int1,int2,int3,...]},...]"

Product
103. getEcNumbersFromProduct()
Input
top

Input Data Type String (for all languages):
email,password
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password)
resultString = client.service.getEcNumbersFromProduct(*parameters)
Output
String containing the different EC Numbers (separated by !) or an array of the different EC Numbers (Python 3) linked to Product, e.g. "EC Number1!EC Number2!EC Number3!EC Number4!..." or "['EC Number1','EC Number2','EC Number3','EC Number4']"

104. getOrganismsFromProduct()
Input
top
Input Data Type String (for all languages):
email,password
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password,"ecNumber*", "product*", "reactionPartners*", "organism*", "ligandStructureId*")
resultString = client.service.getOrganismsFromProduct(*parameters)
Output
String containing the different Organisms (separated by !) or an array of Organisms (Python 3) linked to Product, e.g. "Organism1!Organism2!Organism3!Organism4!..." or "['Organism1','Organism2','Organism3','Organism4']"

105. getProduct(String)
Input
top
Input Data Type String (for all languages):

email,password,"ecNumber*1.1.1.1#organism*Mus musculus"
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password,"ecNumber*", "product*", "reactionPartners*", "organism*", "ligandStructureId*")
resultString = client.service.getProduct(*parameters)
Either the key field ecNumber (e.g. "ecNumber*1.1.1.2"), the key field organism (e.g. "organism*Homo sapiens") or both (e.g. "ecNumber*1.1.1.2#organism*Homo sapiens") have to be specified. If none of these key fields is used as an input parameter, an empty String will be returned by the SOAP query.
In addition, the following optional parameters can be specified: "product", "reactionPartners", "ligandStructureId
Output
String containing the Product entries (entries separated by !, fields separated by #) or an array of Product objects (Python 3), e.g. "ecNumber*string#product*string#reactionPartners*string#organism*string#ligandStructureId*string!
 ecNumber*string#product*string#reactionPartners*string#organism*string#ligandStructureId*string!..."
or
[{'ecNumber':string, 'product':string, 'reactionPartners':string, 'organism':string, 'ligandStructureId':string},{'ecNumber':string, 'product':string, 'reactionPartners':string, 'organism':string, 'ligandStructureId':string},...]"

Purification
106. getEcNumbersFromPurification()
Input
top

Input Data Type String (for all languages):
email,password
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password)
resultString = client.service.getEcNumbersFromPurification(*parameters)
Output
String containing the different EC Numbers (separated by !) or an array of the different EC Numbers (Python 3) linked to Purification, e.g. "EC Number1!EC Number2!EC Number3!EC Number4!..." or "['EC Number1','EC Number2','EC Number3','EC Number4']"

107. getOrganismsFromPurification()
Input
top
Input Data Type String (for all languages):
email,password
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password,"ecNumber*", "commentary*", "organism*", "literature*")
resultString = client.service.getOrganismsFromPurification(*parameters)
Output
String containing the different Organisms (separated by !) or an array of Organisms (Python 3) linked to Purification, e.g. "Organism1!Organism2!Organism3!Organism4!..." or "['Organism1','Organism2','Organism3','Organism4']"

108. getPurification(String)
Input
top
Input Data Type String (for all languages):

email,password,"ecNumber*1.1.1.1#organism*Mus musculus"
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password,"ecNumber*", "commentary*", "organism*", "literature*")
resultString = client.service.getPurification(*parameters)
Either the key field ecNumber (e.g. "ecNumber*1.1.1.2"), the key field organism (e.g. "organism*Homo sapiens") or both (e.g. "ecNumber*1.1.1.2#organism*Homo sapiens") have to be specified. If none of these key fields is used as an input parameter, an empty String will be returned by the SOAP query.
In addition, the following optional parameters can be specified: "commentary", "literature
Output
String containing the Purification entries (entries separated by !, fields separated by #) or an array of Purification objects (Python 3), e.g. "ecNumber*string#commentary*string#organism*string!
 ecNumber*string#commentary*string#organism*string!..."
or
[{'ecNumber':string, 'commentary':string, 'organism':string, 'literature':[int1,int2,int3,...]},{'ecNumber':string, 'commentary':string, 'organism':string, 'literature':[int1,int2,int3,...]},...]"

Reaction
109. getEcNumbersFromReaction()
Input
top

Input Data Type String (for all languages):
email,password
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password)
resultString = client.service.getEcNumbersFromReaction(*parameters)
Output
String containing the different EC Numbers (separated by !) or an array of the different EC Numbers (Python 3) linked to Reaction, e.g. "EC Number1!EC Number2!EC Number3!EC Number4!..." or "['EC Number1','EC Number2','EC Number3','EC Number4']"

110. getOrganismsFromReaction()
Input
top
Input Data Type String (for all languages):
email,password
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password,"ecNumber*", "reaction*", "commentary*", "literature*", "organism*")
resultString = client.service.getOrganismsFromReaction(*parameters)
Output
String containing the different Organisms (separated by !) or an array of Organisms (Python 3) linked to Reaction, e.g. "Organism1!Organism2!Organism3!Organism4!..." or "['Organism1','Organism2','Organism3','Organism4']"

111. getReaction(String)
Input
top
Input Data Type String (for all languages):

email,password,"ecNumber*1.1.1.1#organism*Mus musculus"
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password,"ecNumber*", "reaction*", "commentary*", "literature*", "organism*")
resultString = client.service.getReaction(*parameters)
Either the key field ecNumber (e.g. "ecNumber*1.1.1.2"), the key field organism (e.g. "organism*Homo sapiens") or both (e.g. "ecNumber*1.1.1.2#organism*Homo sapiens") have to be specified. If none of these key fields is used as an input parameter, an empty String will be returned by the SOAP query.
In addition, the following optional parameters can be specified: "reaction", "commentary", "literature
Output
String containing the Reaction entries (entries separated by !, fields separated by #) or an array of Reaction objects (Python 3), e.g. "ecNumber*string#reaction*string#commentary*string#organism*string!
 ecNumber*string#reaction*string#commentary*string#organism*string!..."
or
[{'ecNumber':string, 'reaction':string, 'commentary':string, 'literature':[int1,int2,int3,...], 'organism':string},{'ecNumber':string, 'reaction':string, 'commentary':string, 'literature':[int1,int2,int3,...], 'organism':string},...]"

Reaction Type
112. getEcNumbersFromReactionType()
Input
top

Input Data Type String (for all languages):
email,password
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password)
resultString = client.service.getEcNumbersFromReactionType(*parameters)
Output
String containing the different EC Numbers (separated by !) or an array of the different EC Numbers (Python 3) linked to Reaction Type, e.g. "EC Number1!EC Number2!EC Number3!EC Number4!..." or "['EC Number1','EC Number2','EC Number3','EC Number4']"

113. getOrganismsFromReactionType()
Input
top
Input Data Type String (for all languages):
email,password
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password,"ecNumber*", "reactionType*", "commentary*", "organism*", "literature*")
resultString = client.service.getOrganismsFromReactionType(*parameters)
Output
String containing the different Organisms (separated by !) or an array of Organisms (Python 3) linked to Reaction Type, e.g. "Organism1!Organism2!Organism3!Organism4!..." or "['Organism1','Organism2','Organism3','Organism4']"

114. getReactionType(String)
Input
top
Input Data Type String (for all languages):

email,password,"ecNumber*1.1.1.1#organism*Mus musculus"
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password,"ecNumber*", "reactionType*", "commentary*", "organism*", "literature*")
resultString = client.service.getReactionType(*parameters)
Either the key field ecNumber (e.g. "ecNumber*1.1.1.2"), the key field organism (e.g. "organism*Homo sapiens") or both (e.g. "ecNumber*1.1.1.2#organism*Homo sapiens") have to be specified. If none of these key fields is used as an input parameter, an empty String will be returned by the SOAP query.
In addition, the following optional parameters can be specified: "reactionType", "commentary", "literature
Output
String containing the Reaction Type entries (entries separated by !, fields separated by #) or an array of Reaction Type objects (Python 3), e.g. "ecNumber*string#reactionType*string#commentary*string#organism*string!
 ecNumber*string#reactionType*string#commentary*string#organism*string!..."
or
[{'ecNumber':string, 'reactionType':string, 'commentary':string, 'organism':string, 'literature':[int1,int2,int3,...]},{'ecNumber':string, 'reactionType':string, 'commentary':string, 'organism':string, 'literature':[int1,int2,int3,...]},...]"

Recommended Name
115. getEcNumbersFromRecommendedName()
Input
top

Input Data Type String (for all languages):
email,password
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password)
resultString = client.service.getEcNumbersFromRecommendedName(*parameters)
Output
String containing the different EC Numbers (separated by !) or an array of the different EC Numbers (Python 3) linked to Recommended Name, e.g. "EC Number1!EC Number2!EC Number3!EC Number4!..." or "['EC Number1','EC Number2','EC Number3','EC Number4']"

116. getRecommendedName(String)
Input
top
Input Data Type String (for all languages):

email,password,"ecNumber*1.1.1.1"
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password,"ecNumber*", "recommendedName*", "goNumber*")
resultString = client.service.getRecommendedName(*parameters)
Either the key field ecNumber (e.g. "ecNumber*1.1.1.2"), the key field organism (e.g. "organism*Homo sapiens") or both (e.g. "ecNumber*1.1.1.2#organism*Homo sapiens") have to be specified. If none of these key fields is used as an input parameter, an empty String will be returned by the SOAP query.
In addition, the following optional parameters can be specified: "recommendedName", "goNumber
Output
String containing the Recommended Name entries (entries separated by !, fields separated by #) or an array of Recommended Name objects (Python 3), e.g. "ecNumber*string#recommendedName*string#goNumber*string!
 ecNumber*string#recommendedName*string#goNumber*string!..."
or
[{'ecNumber':string, 'recommendedName':string, 'goNumber':string},{'ecNumber':string, 'recommendedName':string, 'goNumber':string},...]"

Reference
117. getEcNumbersFromReference()
Input
top

Input Data Type String (for all languages):
email,password
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password)
resultString = client.service.getEcNumbersFromReference(*parameters)
Output
String containing the different EC Numbers (separated by !) or an array of the different EC Numbers (Python 3) linked to Reference, e.g. "EC Number1!EC Number2!EC Number3!EC Number4!..." or "['EC Number1','EC Number2','EC Number3','EC Number4']"

118. getOrganismsFromReference()
Input
top
Input Data Type String (for all languages):
email,password
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password,"ecNumber*", "reference*", "authors*", "title*", "journal*", "volume*", "pages*", "year*", "organism*", "commentary*", "pubmedId*", "textmining*")
resultString = client.service.getOrganismsFromReference(*parameters)
Output
String containing the different Organisms (separated by !) or an array of Organisms (Python 3) linked to Reference, e.g. "Organism1!Organism2!Organism3!Organism4!..." or "['Organism1','Organism2','Organism3','Organism4']"

119. getReference(String)
Input
top
Input Data Type String (for all languages):

email,password,"ecNumber*1.1.1.1#organism*Mus musculus"
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password,"ecNumber*", "reference*", "authors*", "title*", "journal*", "volume*", "pages*", "year*", "organism*", "commentary*", "pubmedId*", "textmining*")
resultString = client.service.getReference(*parameters)
Either the key field ecNumber (e.g. "ecNumber*1.1.1.2"), the key field organism (e.g. "organism*Homo sapiens") or both (e.g. "ecNumber*1.1.1.2#organism*Homo sapiens") have to be specified. If none of these key fields is used as an input parameter, an empty String will be returned by the SOAP query.
In addition, the following optional parameters can be specified: "reference", "authors", "title", "journal", "volume", "pages", "year", "commentary", "pubmedId", "textmining
Setting the field/parameter textmining to zero ("textmining*0"), restricts the search to manually annotated entries only.
Output
String containing the Reference entries (entries separated by !, fields separated by #) or an array of Reference objects (Python 3), e.g. "ecNumber*string#reference*string#authors*string#title*string#journal*string#volume*string#pages*string#year*string#organism*string#commentary*string#pubmedId*string#textmining*string!
 ecNumber*string#reference*string#authors*string#title*string#journal*string#volume*string#pages*string#year*string#organism*string#commentary*string#pubmedId*string#textmining*string!..."
or
[{'ecNumber':string, 'reference':string, 'authors':string, 'title':string, 'journal':string, 'volume':string, 'pages':string, 'year':string, 'organism':string, 'commentary':string, 'pubmedId':string, 'textmining':string},{'ecNumber':string, 'reference':string, 'authors':string, 'title':string, 'journal':string, 'volume':string, 'pages':string, 'year':string, 'organism':string, 'commentary':string, 'pubmedId':string, 'textmining':string},...]"

Renatured
120. getEcNumbersFromRenatured()
Input
top

Input Data Type String (for all languages):
email,password
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password)
resultString = client.service.getEcNumbersFromRenatured(*parameters)
Output
String containing the different EC Numbers (separated by !) or an array of the different EC Numbers (Python 3) linked to Renatured, e.g. "EC Number1!EC Number2!EC Number3!EC Number4!..." or "['EC Number1','EC Number2','EC Number3','EC Number4']"

121. getOrganismsFromRenatured()
Input
top
Input Data Type String (for all languages):
email,password
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password,"ecNumber*", "commentary*", "organism*", "literature*")
resultString = client.service.getOrganismsFromRenatured(*parameters)
Output
String containing the different Organisms (separated by !) or an array of Organisms (Python 3) linked to Renatured, e.g. "Organism1!Organism2!Organism3!Organism4!..." or "['Organism1','Organism2','Organism3','Organism4']"

122. getRenatured(String)
Input
top
Input Data Type String (for all languages):

email,password,"ecNumber*1.1.1.1#organism*Mus musculus"
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password,"ecNumber*", "commentary*", "organism*", "literature*")
resultString = client.service.getRenatured(*parameters)
Either the key field ecNumber (e.g. "ecNumber*1.1.1.2"), the key field organism (e.g. "organism*Homo sapiens") or both (e.g. "ecNumber*1.1.1.2#organism*Homo sapiens") have to be specified. If none of these key fields is used as an input parameter, an empty String will be returned by the SOAP query.
In addition, the following optional parameters can be specified: "commentary", "literature
Output
String containing the Renatured entries (entries separated by !, fields separated by #) or an array of Renatured objects (Python 3), e.g. "ecNumber*string#commentary*string#organism*string!
 ecNumber*string#commentary*string#organism*string!..."
or
[{'ecNumber':string, 'commentary':string, 'organism':string, 'literature':[int1,int2,int3,...]},{'ecNumber':string, 'commentary':string, 'organism':string, 'literature':[int1,int2,int3,...]},...]"

Sequence
123. getEcNumbersFromSequence()
Input
top

Input Data Type String (for all languages):
email,password
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password)
resultString = client.service.getEcNumbersFromSequence(*parameters)
Output
String containing the different EC Numbers (separated by !) or an array of the different EC Numbers (Python 3) linked to Sequence, e.g. "EC Number1!EC Number2!EC Number3!EC Number4!..." or "['EC Number1','EC Number2','EC Number3','EC Number4']"

124. getOrganismsFromSequence()
Input
top
Input Data Type String (for all languages):
email,password
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password,"ecNumber*", "sequence*", "noOfAminoAcids*", "firstAccessionCode*", "source*", "id*", "organism*")
resultString = client.service.getOrganismsFromSequence(*parameters)
Output
String containing the different Organisms (separated by !) or an array of Organisms (Python 3) linked to Sequence, e.g. "Organism1!Organism2!Organism3!Organism4!..." or "['Organism1','Organism2','Organism3','Organism4']"

125. getSequence(String)
Input
top
Input Data Type String (for all languages):

email,password,"ecNumber*1.1.1.1#organism*Mus musculus"
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password,"ecNumber*", "sequence*", "noOfAminoAcids*", "firstAccessionCode*", "source*", "id*", "organism*")
resultString = client.service.getSequence(*parameters)
Either the key field ecNumber (e.g. "ecNumber*1.1.1.2"), the key field organism (e.g. "organism*Homo sapiens") or both (e.g. "ecNumber*1.1.1.2#organism*Homo sapiens") have to be specified. If none of these key fields is used as an input parameter, an empty String will be returned by the SOAP query.
In addition, the following optional parameters can be specified: "sequence", "noOfAminoAcids", "firstAccessionCode", "source", "id
Output
String containing the Sequence entries (entries separated by !, fields separated by #) or an array of Sequence objects (Python 3), e.g. "ecNumber*string#sequence*string#noOfAminoAcids*string#firstAccessionCode*string#source*string#id*string#organism*string!
 ecNumber*string#sequence*string#noOfAminoAcids*string#firstAccessionCode*string#source*string#id*string#organism*string!..."
or
[{'ecNumber':string, 'sequence':string, 'noOfAminoAcids':string, 'firstAccessionCode':string, 'source':string, 'id':string, 'organism':string},{'ecNumber':string, 'sequence':string, 'noOfAminoAcids':string, 'firstAccessionCode':string, 'source':string, 'id':string, 'organism':string},...]"

Source Tissue
126. getEcNumbersFromSourceTissue()
Input
top

Input Data Type String (for all languages):
email,password
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password)
resultString = client.service.getEcNumbersFromSourceTissue(*parameters)
Output
String containing the different EC Numbers (separated by !) or an array of the different EC Numbers (Python 3) linked to Source Tissue, e.g. "EC Number1!EC Number2!EC Number3!EC Number4!..." or "['EC Number1','EC Number2','EC Number3','EC Number4']"

127. getOrganismsFromSourceTissue()
Input
top
Input Data Type String (for all languages):
email,password
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password,"ecNumber*", "sourceTissue*", "commentary*", "organism*", "literature*", "textmining*")
resultString = client.service.getOrganismsFromSourceTissue(*parameters)
Output
String containing the different Organisms (separated by !) or an array of Organisms (Python 3) linked to Source Tissue, e.g. "Organism1!Organism2!Organism3!Organism4!..." or "['Organism1','Organism2','Organism3','Organism4']"

128. getSourceTissue(String)
Input
top
Input Data Type String (for all languages):

email,password,"ecNumber*1.1.1.1#organism*Mus musculus"
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password,"ecNumber*", "sourceTissue*", "commentary*", "organism*", "literature*", "textmining*")
resultString = client.service.getSourceTissue(*parameters)
Either the key field ecNumber (e.g. "ecNumber*1.1.1.2"), the key field organism (e.g. "organism*Homo sapiens") or both (e.g. "ecNumber*1.1.1.2#organism*Homo sapiens") have to be specified. If none of these key fields is used as an input parameter, an empty String will be returned by the SOAP query.
In addition, the following optional parameters can be specified: "sourceTissue", "commentary", "literature", "textmining
Setting the field/parameter textmining to zero ("textmining*0"), restricts the search to manually annotated entries only.
Output
String containing the Source Tissue entries (entries separated by !, fields separated by #) or an array of Source Tissue objects (Python 3), e.g. "ecNumber*string#sourceTissue*string#commentary*string#organism*string#textmining*string!
 ecNumber*string#sourceTissue*string#commentary*string#organism*string#textmining*string!..."
or
[{'ecNumber':string, 'sourceTissue':string, 'commentary':string, 'organism':string, 'literature':[int1,int2,int3,...], 'textmining':string},{'ecNumber':string, 'sourceTissue':string, 'commentary':string, 'organism':string, 'literature':[int1,int2,int3,...], 'textmining':string},...]"

Specific Activity
129. getEcNumbersFromSpecificActivity()
Input
top

Input Data Type String (for all languages):
email,password
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password)
resultString = client.service.getEcNumbersFromSpecificActivity(*parameters)
Output
String containing the different EC Numbers (separated by !) or an array of the different EC Numbers (Python 3) linked to Specific Activity, e.g. "EC Number1!EC Number2!EC Number3!EC Number4!..." or "['EC Number1','EC Number2','EC Number3','EC Number4']"

130. getOrganismsFromSpecificActivity()
Input
top
Input Data Type String (for all languages):
email,password
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password,"ecNumber*", "specificActivity*", "specificActivityMaximum*", "commentary*", "organism*", "literature*")
resultString = client.service.getOrganismsFromSpecificActivity(*parameters)
Output
String containing the different Organisms (separated by !) or an array of Organisms (Python 3) linked to Specific Activity, e.g. "Organism1!Organism2!Organism3!Organism4!..." or "['Organism1','Organism2','Organism3','Organism4']"

131. getSpecificActivity(String)
Input
top
Input Data Type String (for all languages):

email,password,"ecNumber*1.1.1.1#organism*Mus musculus"
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password,"ecNumber*", "specificActivity*", "specificActivityMaximum*", "commentary*", "organism*", "literature*")
resultString = client.service.getSpecificActivity(*parameters)
Either the key field ecNumber (e.g. "ecNumber*1.1.1.2"), the key field organism (e.g. "organism*Homo sapiens") or both (e.g. "ecNumber*1.1.1.2#organism*Homo sapiens") have to be specified. If none of these key fields is used as an input parameter, an empty String will be returned by the SOAP query.
In addition, the following optional parameters can be specified: "specificActivity", "specificActivityMaximum", "commentary", "literature
Output
String containing the Specific Activity entries (entries separated by !, fields separated by #) or an array of Specific Activity objects (Python 3), e.g. "ecNumber*string#specificActivity*string#specificActivityMaximum*string#commentary*string#organism*string!
 ecNumber*string#specificActivity*string#specificActivityMaximum*string#commentary*string#organism*string!..."
or
[{'ecNumber':string, 'specificActivity':string, 'specificActivityMaximum':string, 'commentary':string, 'organism':string, 'literature':[int1,int2,int3,...]},{'ecNumber':string, 'specificActivity':string, 'specificActivityMaximum':string, 'commentary':string, 'organism':string, 'literature':[int1,int2,int3,...]},...]"

Storage Stability
132. getEcNumbersFromStorageStability()
Input
top

Input Data Type String (for all languages):
email,password
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password)
resultString = client.service.getEcNumbersFromStorageStability(*parameters)
Output
String containing the different EC Numbers (separated by !) or an array of the different EC Numbers (Python 3) linked to Storage Stability, e.g. "EC Number1!EC Number2!EC Number3!EC Number4!..." or "['EC Number1','EC Number2','EC Number3','EC Number4']"

133. getOrganismsFromStorageStability()
Input
top
Input Data Type String (for all languages):
email,password
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password,"ecNumber*", "storageStability*", "organism*", "literature*")
resultString = client.service.getOrganismsFromStorageStability(*parameters)
Output
String containing the different Organisms (separated by !) or an array of Organisms (Python 3) linked to Storage Stability, e.g. "Organism1!Organism2!Organism3!Organism4!..." or "['Organism1','Organism2','Organism3','Organism4']"

134. getStorageStability(String)
Input
top
Input Data Type String (for all languages):

email,password,"ecNumber*1.1.1.1#organism*Mus musculus"
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password,"ecNumber*", "storageStability*", "organism*", "literature*")
resultString = client.service.getStorageStability(*parameters)
Either the key field ecNumber (e.g. "ecNumber*1.1.1.2"), the key field organism (e.g. "organism*Homo sapiens") or both (e.g. "ecNumber*1.1.1.2#organism*Homo sapiens") have to be specified. If none of these key fields is used as an input parameter, an empty String will be returned by the SOAP query.
In addition, the following optional parameters can be specified: "storageStability", "literature
Output
String containing the Storage Stability entries (entries separated by !, fields separated by #) or an array of Storage Stability objects (Python 3), e.g. "ecNumber*string#storageStability*string#organism*string!
 ecNumber*string#storageStability*string#organism*string!..."
or
[{'ecNumber':string, 'storageStability':string, 'organism':string, 'literature':[int1,int2,int3,...]},{'ecNumber':string, 'storageStability':string, 'organism':string, 'literature':[int1,int2,int3,...]},...]"

Substrate
135. getEcNumbersFromSubstrate()
Input
top

Input Data Type String (for all languages):
email,password
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password)
resultString = client.service.getEcNumbersFromSubstrate(*parameters)
Output
String containing the different EC Numbers (separated by !) or an array of the different EC Numbers (Python 3) linked to Substrate, e.g. "EC Number1!EC Number2!EC Number3!EC Number4!..." or "['EC Number1','EC Number2','EC Number3','EC Number4']"

136. getOrganismsFromSubstrate()
Input
top
Input Data Type String (for all languages):
email,password
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password,"ecNumber*", "substrate*", "reactionPartners*", "organism*", "ligandStructureId*")
resultString = client.service.getOrganismsFromSubstrate(*parameters)
Output
String containing the different Organisms (separated by !) or an array of Organisms (Python 3) linked to Substrate, e.g. "Organism1!Organism2!Organism3!Organism4!..." or "['Organism1','Organism2','Organism3','Organism4']"

137. getSubstrate(String)
Input
top
Input Data Type String (for all languages):

email,password,"ecNumber*1.1.1.1#organism*Mus musculus"
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password,"ecNumber*", "substrate*", "reactionPartners*", "organism*", "ligandStructureId*")
resultString = client.service.getSubstrate(*parameters)
Either the key field ecNumber (e.g. "ecNumber*1.1.1.2"), the key field organism (e.g. "organism*Homo sapiens") or both (e.g. "ecNumber*1.1.1.2#organism*Homo sapiens") have to be specified. If none of these key fields is used as an input parameter, an empty String will be returned by the SOAP query.
In addition, the following optional parameters can be specified: "substrate", "reactionPartners", "ligandStructureId
Output
String containing the Substrate entries (entries separated by !, fields separated by #) or an array of Substrate objects (Python 3), e.g. "ecNumber*string#substrate*string#reactionPartners*string#organism*string#ligandStructureId*string!
 ecNumber*string#substrate*string#reactionPartners*string#organism*string#ligandStructureId*string!..."
or
[{'ecNumber':string, 'substrate':string, 'reactionPartners':string, 'organism':string, 'ligandStructureId':string},{'ecNumber':string, 'substrate':string, 'reactionPartners':string, 'organism':string, 'ligandStructureId':string},...]"

Substrates Products
138. getEcNumbersFromSubstratesProducts()
Input
top

Input Data Type String (for all languages):
email,password
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password)
resultString = client.service.getEcNumbersFromSubstratesProducts(*parameters)
Output
String containing the different EC Numbers (separated by !) or an array of the different EC Numbers (Python 3) linked to Substrates Products, e.g. "EC Number1!EC Number2!EC Number3!EC Number4!..." or "['EC Number1','EC Number2','EC Number3','EC Number4']"

139. getSubstratesProducts(String)
Input
top
Input Data Type String (for all languages):

email,password,"ecNumber*1.1.1.1"
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password,"ecNumber*", "substrates*", "commentarySubstrates*", "literatureSubstrates*", "organismSubstrates*", "products*", "commentaryProducts*", "literatureProducts*", "organismProducts*", "reversibility*")
resultString = client.service.getSubstratesProducts(*parameters)
Either the key field ecNumber (e.g. "ecNumber*1.1.1.2"), the key field organism (e.g. "organism*Homo sapiens") or both (e.g. "ecNumber*1.1.1.2#organism*Homo sapiens") have to be specified. If none of these key fields is used as an input parameter, an empty String will be returned by the SOAP query.
In addition, the following optional parameters can be specified: "substrates", "commentarySubstrates", "literatureSubstrates", "organismSubstrates", "products", "commentaryProducts", "literatureProducts", "organismProducts", "reversibility
Output
String containing the Substrates Products entries (entries separated by !, fields separated by #) or an array of Substrates Products objects (Python 3), e.g. "ecNumber*string#substrates*string#commentarySubstrates*string#literatureSubstrates*string#organismSubstrates*string#products*string#commentaryProducts*string#literatureProducts*string#organismProducts*string#reversibility*string!
 ecNumber*string#substrates*string#commentarySubstrates*string#literatureSubstrates*string#organismSubstrates*string#products*string#commentaryProducts*string#literatureProducts*string#organismProducts*string#reversibility*string!..."
or
[{'ecNumber':string, 'substrates':string, 'commentarySubstrates':string, 'literatureSubstrates':string, 'organismSubstrates':string, 'products':string, 'commentaryProducts':string, 'literatureProducts':string, 'organismProducts':string, 'reversibility':string},{'ecNumber':string, 'substrates':string, 'commentarySubstrates':string, 'literatureSubstrates':string, 'organismSubstrates':string, 'products':string, 'commentaryProducts':string, 'literatureProducts':string, 'organismProducts':string, 'reversibility':string},...]"

Subunits
140. getEcNumbersFromSubunits()
Input
top

Input Data Type String (for all languages):
email,password
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password)
resultString = client.service.getEcNumbersFromSubunits(*parameters)
Output
String containing the different EC Numbers (separated by !) or an array of the different EC Numbers (Python 3) linked to Subunits, e.g. "EC Number1!EC Number2!EC Number3!EC Number4!..." or "['EC Number1','EC Number2','EC Number3','EC Number4']"

141. getOrganismsFromSubunits()
Input
top
Input Data Type String (for all languages):
email,password
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password,"ecNumber*", "subunits*", "commentary*", "organism*", "literature*")
resultString = client.service.getOrganismsFromSubunits(*parameters)
Output
String containing the different Organisms (separated by !) or an array of Organisms (Python 3) linked to Subunits, e.g. "Organism1!Organism2!Organism3!Organism4!..." or "['Organism1','Organism2','Organism3','Organism4']"

142. getSubunits(String)
Input
top
Input Data Type String (for all languages):

email,password,"ecNumber*1.1.1.1#organism*Mus musculus"
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password,"ecNumber*", "subunits*", "commentary*", "organism*", "literature*")
resultString = client.service.getSubunits(*parameters)
Either the key field ecNumber (e.g. "ecNumber*1.1.1.2"), the key field organism (e.g. "organism*Homo sapiens") or both (e.g. "ecNumber*1.1.1.2#organism*Homo sapiens") have to be specified. If none of these key fields is used as an input parameter, an empty String will be returned by the SOAP query.
In addition, the following optional parameters can be specified: "subunits", "commentary", "literature
Output
String containing the Subunits entries (entries separated by !, fields separated by #) or an array of Subunits objects (Python 3), e.g. "ecNumber*string#subunits*string#commentary*string#organism*string!
 ecNumber*string#subunits*string#commentary*string#organism*string!..."
or
[{'ecNumber':string, 'subunits':string, 'commentary':string, 'organism':string, 'literature':[int1,int2,int3,...]},{'ecNumber':string, 'subunits':string, 'commentary':string, 'organism':string, 'literature':[int1,int2,int3,...]},...]"

Synonyms
143. getEcNumbersFromSynonyms()
Input
top

Input Data Type String (for all languages):
email,password
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password)
resultString = client.service.getEcNumbersFromSynonyms(*parameters)
Output
String containing the different EC Numbers (separated by !) or an array of the different EC Numbers (Python 3) linked to Synonyms, e.g. "EC Number1!EC Number2!EC Number3!EC Number4!..." or "['EC Number1','EC Number2','EC Number3','EC Number4']"

144. getOrganismsFromSynonyms()
Input
top
Input Data Type String (for all languages):
email,password
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password,"ecNumber*", "synonyms*", "commentary*", "organism*", "literature*")
resultString = client.service.getOrganismsFromSynonyms(*parameters)
Output
String containing the different Organisms (separated by !) or an array of Organisms (Python 3) linked to Synonyms, e.g. "Organism1!Organism2!Organism3!Organism4!..." or "['Organism1','Organism2','Organism3','Organism4']"

145. getSynonyms(String)
Input
top
Input Data Type String (for all languages):

email,password,"ecNumber*1.1.1.1#organism*Mus musculus"
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password,"ecNumber*", "synonyms*", "commentary*", "organism*", "literature*")
resultString = client.service.getSynonyms(*parameters)
Either the key field ecNumber (e.g. "ecNumber*1.1.1.2"), the key field organism (e.g. "organism*Homo sapiens") or both (e.g. "ecNumber*1.1.1.2#organism*Homo sapiens") have to be specified. If none of these key fields is used as an input parameter, an empty String will be returned by the SOAP query.
In addition, the following optional parameters can be specified: "synonyms", "commentary", "literature
Output
String containing the Synonyms entries (entries separated by !, fields separated by #) or an array of Synonyms objects (Python 3), e.g. "ecNumber*string#synonyms*string#commentary*string#organism*string!
 ecNumber*string#synonyms*string#commentary*string#organism*string!..."
or
[{'ecNumber':string, 'synonyms':string, 'commentary':string, 'organism':string, 'literature':[int1,int2,int3,...]},{'ecNumber':string, 'synonyms':string, 'commentary':string, 'organism':string, 'literature':[int1,int2,int3,...]},...]"

Systematic Name
146. getEcNumbersFromSystematicName()
Input
top

Input Data Type String (for all languages):
email,password
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password)
resultString = client.service.getEcNumbersFromSystematicName(*parameters)
Output
String containing the different EC Numbers (separated by !) or an array of the different EC Numbers (Python 3) linked to Systematic Name, e.g. "EC Number1!EC Number2!EC Number3!EC Number4!..." or "['EC Number1','EC Number2','EC Number3','EC Number4']"

147. getSystematicName(String)
Input
top
Input Data Type String (for all languages):

email,password,"ecNumber*1.1.1.1"
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password,"ecNumber*", "systematicName*")
resultString = client.service.getSystematicName(*parameters)
Either the key field ecNumber (e.g. "ecNumber*1.1.1.2"), the key field organism (e.g. "organism*Homo sapiens") or both (e.g. "ecNumber*1.1.1.2#organism*Homo sapiens") have to be specified. If none of these key fields is used as an input parameter, an empty String will be returned by the SOAP query.
In addition, the following optional parameters can be specified: "systematicName
Output
String containing the Systematic Name entries (entries separated by !, fields separated by #) or an array of Systematic Name objects (Python 3), e.g. "ecNumber*string#systematicName*string!
 ecNumber*string#systematicName*string!..."
or
[{'ecNumber':string, 'systematicName':string},{'ecNumber':string, 'systematicName':string},...]"

Temperature Optimum
148. getEcNumbersFromTemperatureOptimum()
Input
top

Input Data Type String (for all languages):
email,password
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password)
resultString = client.service.getEcNumbersFromTemperatureOptimum(*parameters)
Output
String containing the different EC Numbers (separated by !) or an array of the different EC Numbers (Python 3) linked to Temperature Optimum, e.g. "EC Number1!EC Number2!EC Number3!EC Number4!..." or "['EC Number1','EC Number2','EC Number3','EC Number4']"

149. getOrganismsFromTemperatureOptimum()
Input
top
Input Data Type String (for all languages):
email,password
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password,"ecNumber*", "temperatureOptimum*", "temperatureOptimumMaximum*", "commentary*", "organism*", "literature*")
resultString = client.service.getOrganismsFromTemperatureOptimum(*parameters)
Output
String containing the different Organisms (separated by !) or an array of Organisms (Python 3) linked to Temperature Optimum, e.g. "Organism1!Organism2!Organism3!Organism4!..." or "['Organism1','Organism2','Organism3','Organism4']"

150. getTemperatureOptimum(String)
Input
top
Input Data Type String (for all languages):

email,password,"ecNumber*1.1.1.1#organism*Mus musculus"
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password,"ecNumber*", "temperatureOptimum*", "temperatureOptimumMaximum*", "commentary*", "organism*", "literature*")
resultString = client.service.getTemperatureOptimum(*parameters)
Either the key field ecNumber (e.g. "ecNumber*1.1.1.2"), the key field organism (e.g. "organism*Homo sapiens") or both (e.g. "ecNumber*1.1.1.2#organism*Homo sapiens") have to be specified. If none of these key fields is used as an input parameter, an empty String will be returned by the SOAP query.
In addition, the following optional parameters can be specified: "temperatureOptimum", "temperatureOptimumMaximum", "commentary", "literature
Output
String containing the Temperature Optimum entries (entries separated by !, fields separated by #) or an array of Temperature Optimum objects (Python 3), e.g. "ecNumber*string#temperatureOptimum*string#temperatureOptimumMaximum*string#commentary*string#organism*string!
 ecNumber*string#temperatureOptimum*string#temperatureOptimumMaximum*string#commentary*string#organism*string!..."
or
[{'ecNumber':string, 'temperatureOptimum':string, 'temperatureOptimumMaximum':string, 'commentary':string, 'organism':string, 'literature':[int1,int2,int3,...]},{'ecNumber':string, 'temperatureOptimum':string, 'temperatureOptimumMaximum':string, 'commentary':string, 'organism':string, 'literature':[int1,int2,int3,...]},...]"

Temperature Range
151. getEcNumbersFromTemperatureRange()
Input
top

Input Data Type String (for all languages):
email,password
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password)
resultString = client.service.getEcNumbersFromTemperatureRange(*parameters)
Output
String containing the different EC Numbers (separated by !) or an array of the different EC Numbers (Python 3) linked to Temperature Range, e.g. "EC Number1!EC Number2!EC Number3!EC Number4!..." or "['EC Number1','EC Number2','EC Number3','EC Number4']"

152. getOrganismsFromTemperatureRange()
Input
top
Input Data Type String (for all languages):
email,password
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password,"ecNumber*", "temperatureRange*", "temperatureRangeMaximum*", "commentary*", "organism*", "literature*")
resultString = client.service.getOrganismsFromTemperatureRange(*parameters)
Output
String containing the different Organisms (separated by !) or an array of Organisms (Python 3) linked to Temperature Range, e.g. "Organism1!Organism2!Organism3!Organism4!..." or "['Organism1','Organism2','Organism3','Organism4']"

153. getTemperatureRange(String)
Input
top
Input Data Type String (for all languages):

email,password,"ecNumber*1.1.1.1#organism*Mus musculus"
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password,"ecNumber*", "temperatureRange*", "temperatureRangeMaximum*", "commentary*", "organism*", "literature*")
resultString = client.service.getTemperatureRange(*parameters)
Either the key field ecNumber (e.g. "ecNumber*1.1.1.2"), the key field organism (e.g. "organism*Homo sapiens") or both (e.g. "ecNumber*1.1.1.2#organism*Homo sapiens") have to be specified. If none of these key fields is used as an input parameter, an empty String will be returned by the SOAP query.
In addition, the following optional parameters can be specified: "temperatureRange", "temperatureRangeMaximum", "commentary", "literature
Output
String containing the Temperature Range entries (entries separated by !, fields separated by #) or an array of Temperature Range objects (Python 3), e.g. "ecNumber*string#temperatureRange*string#temperatureRangeMaximum*string#commentary*string#organism*string!
 ecNumber*string#temperatureRange*string#temperatureRangeMaximum*string#commentary*string#organism*string!..."
or
[{'ecNumber':string, 'temperatureRange':string, 'temperatureRangeMaximum':string, 'commentary':string, 'organism':string, 'literature':[int1,int2,int3,...]},{'ecNumber':string, 'temperatureRange':string, 'temperatureRangeMaximum':string, 'commentary':string, 'organism':string, 'literature':[int1,int2,int3,...]},...]"

Temperature Stability
154. getEcNumbersFromTemperatureStability()
Input
top

Input Data Type String (for all languages):
email,password
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password)
resultString = client.service.getEcNumbersFromTemperatureStability(*parameters)
Output
String containing the different EC Numbers (separated by !) or an array of the different EC Numbers (Python 3) linked to Temperature Stability, e.g. "EC Number1!EC Number2!EC Number3!EC Number4!..." or "['EC Number1','EC Number2','EC Number3','EC Number4']"

155. getOrganismsFromTemperatureStability()
Input
top
Input Data Type String (for all languages):
email,password
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password,"ecNumber*", "temperatureStability*", "temperatureStabilityMaximum*", "commentary*", "organism*", "literature*")
resultString = client.service.getOrganismsFromTemperatureStability(*parameters)
Output
String containing the different Organisms (separated by !) or an array of Organisms (Python 3) linked to Temperature Stability, e.g. "Organism1!Organism2!Organism3!Organism4!..." or "['Organism1','Organism2','Organism3','Organism4']"

156. getTemperatureStability(String)
Input
top
Input Data Type String (for all languages):

email,password,"ecNumber*1.1.1.1#organism*Mus musculus"
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password,"ecNumber*", "temperatureStability*", "temperatureStabilityMaximum*", "commentary*", "organism*", "literature*")
resultString = client.service.getTemperatureStability(*parameters)
Either the key field ecNumber (e.g. "ecNumber*1.1.1.2"), the key field organism (e.g. "organism*Homo sapiens") or both (e.g. "ecNumber*1.1.1.2#organism*Homo sapiens") have to be specified. If none of these key fields is used as an input parameter, an empty String will be returned by the SOAP query.
In addition, the following optional parameters can be specified: "temperatureStability", "temperatureStabilityMaximum", "commentary", "literature
Output
String containing the Temperature Stability entries (entries separated by !, fields separated by #) or an array of Temperature Stability objects (Python 3), e.g. "ecNumber*string#temperatureStability*string#temperatureStabilityMaximum*string#commentary*string#organism*string!
 ecNumber*string#temperatureStability*string#temperatureStabilityMaximum*string#commentary*string#organism*string!..."
or
[{'ecNumber':string, 'temperatureStability':string, 'temperatureStabilityMaximum':string, 'commentary':string, 'organism':string, 'literature':[int1,int2,int3,...]},{'ecNumber':string, 'temperatureStability':string, 'temperatureStabilityMaximum':string, 'commentary':string, 'organism':string, 'literature':[int1,int2,int3,...]},...]"

Turnover Number
157. getEcNumbersFromTurnoverNumber()
Input
top

Input Data Type String (for all languages):
email,password
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password)
resultString = client.service.getEcNumbersFromTurnoverNumber(*parameters)
Output
String containing the different EC Numbers (separated by !) or an array of the different EC Numbers (Python 3) linked to Turnover Number, e.g. "EC Number1!EC Number2!EC Number3!EC Number4!..." or "['EC Number1','EC Number2','EC Number3','EC Number4']"

158. getOrganismsFromTurnoverNumber()
Input
top
Input Data Type String (for all languages):
email,password
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password,"ecNumber*", "turnoverNumber*", "turnoverNumberMaximum*", "substrate*", "commentary*", "organism*", "ligandStructureId*", "literature*")
resultString = client.service.getOrganismsFromTurnoverNumber(*parameters)
Output
String containing the different Organisms (separated by !) or an array of Organisms (Python 3) linked to Turnover Number, e.g. "Organism1!Organism2!Organism3!Organism4!..." or "['Organism1','Organism2','Organism3','Organism4']"

159. getTurnoverNumber(String)
Input
top
Input Data Type String (for all languages):

email,password,"ecNumber*1.1.1.1#organism*Mus musculus"
Python 3 example code snippet:
parameters = ("j.doe@example.edu",password,"ecNumber*", "turnoverNumber*", "turnoverNumberMaximum*", "substrate*", "commentary*", "organism*", "ligandStructureId*", "literature*")
resultString = client.service.getTurnoverNumber(*parameters)
Either the key field ecNumber (e.g. "ecNumber*1.1.1.2"), the key field organism (e.g. "organism*Homo sapiens") or both (e.g. "ecNumber*1.1.1.2#organism*Homo sapiens") have to be specified. If none of these key fields is used as an input parameter, an empty String will be returned by the SOAP query.
In addition, the following optional parameters can be specified: "turnoverNumber", "turnoverNumberMaximum", "substrate", "commentary", "ligandStructureId", "literature
Output
String containing the Turnover Number entries (entries separated by !, fields separated by #) or an array of Turnover Number objects (Python 3), e.g. "ecNumber*string#turnoverNumber*string#turnoverNumberMaximum*string#substrate*string#commentary*string#organism*string#ligandStructureId*string!
 ecNumber*string#turnoverNumber*string#turnoverNumberMaximum*string#substrate*string#commentary*string#organism*string#ligandStructureId*string!..."
or
[{'ecNumber':string, 'turnoverNumber':string, 'turnoverNumberMaximum':string, 'substrate':string, 'commentary':string, 'organism':string, 'ligandStructureId':string, 'literature':[int1,int2,int3,...]},{'ecNumber':string, 'turnoverNumber':string, 'turnoverNumberMaximum':string, 'substrate':string, 'commentary':string, 'organism':string, 'ligandStructureId':string, 'literature':[int1,int2,int3,...]},...]"