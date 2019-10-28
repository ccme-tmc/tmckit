#!/usr/bin/env python
#This library define XML schema (xsd) based XML I/O for generic python routine
#By using xsd, it is possible to convert between xml file and Python class object, or list object with definition only thus standarize the I/O.
#The XML schema follows standard, and some parameters have special meaning here.
#Schema format:
#In XSD, there must be one ( and only first is used ) xs:element in direct child of root in schema followed by a lot of definitions. The xs:element represent the created object.
#In XML, there must be one ( and only one ) node in direct children of the root, corresponding to the object. Note, the name of the top node is not used as there must be only one. However it is suggest to give it a meaningful name; when serilization it is set to the class name
#name: name is the type name, must correspoding to a class name in python code. If start as ArrayOf, this object will be a list object; no class definition is needed here.
#minOccurs/maxOccurs : Due to the nature of class/list, we do not present methods to deal with variable-length objects, so minOccurs must equals to maxOccurs. There are two exceptions:
#1. the last element in list type; the maxOccurs can be "unbounded", means everything follows belong to this.
#2. any element in class type: minOccurs can be 0, which means the object can be neglected in XML
#strings are converted to ANSI string ( they are Unicode in XML)
#
#Speical type: ArrayOfX
#This means thie type is stored in list instead of class
#
#Special type: TableOfX
#This means the whole list is stored in a space/tab+endline-seperated string
#It must has only exactly one element type is TableRowOfX
#It must correspond to a list type

#Special type: TableRowOfX
#This means the whole class/list is stored in a space/tab-seperated string
#It can be either list or class(which actual class name is X). If it is a list, its name must be ArrayOfTableRowOfX

#Warning: Use TableOf will cause XML cannot be validated by XSD! Please be cautious when using it.

from __future__ import print_function
import xml.dom.minidom as minidom
import os

def f_xml_GetSchemaFilePath(stFileName):
    '''
    return full schema file path from file name
    '''
    return os.path.join(os.path.dirname(__file__),"schemas/"+stFileName) 

def f_xml_GetChildElement(node):
    '''
    Get direct child element of one node ( to remove all text node )
    '''
    return [x for x in node.childNodes if x.nodeType == x.ELEMENT_NODE] 

def f_xml_GetDirectChildElementByNodeName(node,stName):
    '''
    Get the direct child element with specific node name
    '''
    return [x for x in node.childNodes if x.nodeType == x.ELEMENT_NODE and x.nodeName == stName] 


def f_xml_IsBaseType(stType):
    '''
    Determine whether a string is an allowed base type.
    Any attribute must be a base type.
    '''
    if ( stType[2] == ":"):
        return True
    else:
        return False

def f_xml_ReadBaseTypeString(value,stType):
    '''
    Convert a string to a base type in XML schema, like "10" to 10, "true" to True.
    :param stType: type name, can be xs:string or string
    '''
    obj = None

    if ( stType[2] == ":" ):
        stType2 = stType[3:]
    else:
        stType2 = stType

    if ( stType2 == "string"):
        obj = str(value)
    elif ( stType2 == "int" or stType2 == "integer"):
        obj = int(value)
    elif ( stType2 == "decimal" ):
        obj = float(value)
    elif ( stType2 == "boolean"):
        if ( value.lower() == "true"):
            obj = True
        elif ( value.lower() == "false"):
            obj = False
        else:
            raise ValueError("Boolean type must has value 'true' or 'false'")
    else:
        raise ValueError("Unknown type: "+stType2)
    return obj

def f_xml_ReadBaseTypeNode(node_xml,stType):
    '''
    Convert a xml node of a base type to object
    :return: object if can be convert, None if nothing in the node_xml
    '''
    if ( len(node_xml.childNodes) == 1): #Text node content
        value = node_xml.childNodes[0].data
        return f_xml_ReadBaseTypeString(value,stType)
    else:
#For string return an empty string, else return None
        if ( stType == "xs:string"):
            return ""
        else:
            return None

 
def f_xml_TrimAllText(stXML):
    '''
    Reformat 
    <a>
    xxx
    </a>
    to 
    <a>xxx</a>
    by replacing all space between ">" and non-"<" or non-">" and "<"
    space, tab and line end are treated as space
    '''
    #print("Reformatt")
    nStart = 0
    listDel=[] #Store all space start and end index for later delete
    while (True):
        n1 = stXML.find(">",nStart)
        n2 = stXML.find("<",n1)
        if ( n1 == -1 or n2 == -1):
            break
        #If it is already a multi-line text, then do not remove end of line 
        n = stXML.count("\n",n1,n2)
        st_remove = " \t\n"
        #Remove start. Do not remove first and final end-of-line if it is multi-line text.
        for i in range(n1+1,n2-1): #Non-Tag Text between > and <
            if ( not stXML[i] in st_remove ):
                if ( i==n1+1 ): #Nothing to remove
                    break
                if ( n == 2):
                    listDel.append([n1+1,i])
                else:
                    listDel.append([n1+2,i])
                break

        for i in range(n2-1,n1+1,-1):
            if (not stXML[i] in st_remove ):
                if (n == 2 ):
                #remove final space for single-line text between two tag
                    if (i != n2-1 ):#Nothing to remove
                        listDel.append([i+1,n2])
                #remove totally empty line
                elif (n > 2):#Remove totally empty line
                    num_line = stXML.count("\n",i,n2-1)
                    if (num_line >1):
                        listDel.append([stXML.index("\n",i)+1,stXML.rindex("\n",i,n2-1)+1])
                break
                        
        nStart = n2

    listDel.reverse()
    for l2 in listDel:
        stXML = stXML[:l2[0]] + stXML[l2[1]:]

    return stXML



class XMLSerilizer():
    '''
    XML schema-based engine to serilize object to xml
    Efficiency is rather low, however we will not parse too large file here.
    '''

    def __init__(self,varSet=[],stSchemeFileName=""):
        '''
        Initialize
        :param varSet: the array contains all class defition may used in xml parsing procedure. Normally, just pass locals() or globals() to it.
        :param stSchemeFileName: the file name of XML scheme
        '''
        self.classDef =varSet
        self.schema = None
        if ( stSchemeFileName != None):
            self.LoadXSD(stSchemeFileName)

    def LoadXSD(self,stFile):
        '''
        Read XML schema for later converting bewteen object and XML file 
        Multiple schema is allowed and they will be combined.
        '''
        if ( not os.path.exists(stFile)):#Search order : current folder, default folder ( tmckit/schemas )
            stFile = f_xml_GetSchemaFilePath(stFile) 
            if ( os.path.exists(f_xml_GetSchemaFilePath(stFile))):
                raise ValueError("schema %s not found" % stFile)
        
        if ( self.schema == None):
            self.schema = f_xml_GetChildElement(minidom.parse(stFile))[0] #Only first schema is used here
        else:
            schema2 = f_xml_GetChildElement(minidom.parse(stFile))[0] #Only first schema is used here
            for childNode in schema2:
                self.schema.append(childNode)

    def __IsArrayElement__(self,node_xsd):
        '''
        Determine whether a element is an array-type element and must be stored in an array
        '''
        if ( node_xsd.hasAttribute("maxOccurs")):
            st = node_xsd.getAttribute("maxOccurs")
            if ( st == "unbounded"):
                return True
            else:
                n = int(st)
                if ( n >= 2):# 1 can also be an array, but 2 must be
                    return True
        return False


    def __IsArray__(self,node_xsd):
        '''
        Determine whether a element is an array type
        A Array ( List) type must be constructed as follows:
        Typename : ArrayOf* ( Any Type named ArrayOf* must follow next conditino )
        Content: every elements must have maxOccurs attribute, like 1, 2 or  "unbounded"
        Array instance definition ( element of other type ), array definition (or array element definition -- this does not work) can be passed to this function 
        '''
        nodeName =  node_xsd.getAttribute("name")
        #print(node_xsd.nodeName,nodeName)
        if ( node_xsd.nodeName == "xs:complexType"): #Array definition
            if ( nodeName[:7] == "ArrayOf" or nodeName[:7] == "TableOf" ):
                childNodes = node_xsd.getElementsByTagName("xs:element")
#Below line is removed, as ArrayOf can be used as storage method indicator ( store a class linear in in list )
                #if ( len(childNodes) > 1):
                #    raise ValueError,"Array type definition can only have one ielement: " + nodeName

#               if ( not self.__IsArrayElement__(childNodes[0])):
                for childNode in childNodes:
                    if ( not childNode.hasAttribute("maxOccurs")):
                        raise ValueError("Element in Array type definition must have maxOccurs attribute: " + nodeName)
                return True
            else:
                return False
        elif ( node_xsd.nodeName =="xs:element"): #Array element
            return self.__IsArray__(self.__FindTypeSchema__(node_xsd.getAttribute("type")))


    def __IsTable__(self,node_xsd):
        '''
        Determine whether an element is an Table or TableRow type
        Only by name
        '''
        nodeName =  node_xsd.getAttribute("name")
        if ( node_xsd.nodeName == "xs:complexType"): #Array definition
            if ( nodeName[:7] == "TableOf" or "TableRowOf" in nodeName ):
                return True
            else:
                #print("%s is not Table" % node_xsd.toxml())
                return False
        elif ( node_xsd.nodeName =="xs:element"): #Array element
            return self.__IsTable__(self.__FindTypeSchema__(node_xsd.getAttribute("type")))


    def __FindTypeSchema__(self,stType):
        '''
        Return the node represent complex type in schema
        '''
        node_xsd = None
        for childNode in f_xml_GetChildElement(self.schema):
            if ( childNode.nodeName == "xs:complexType"):
                if ( childNode.getAttribute("name") == stType):
                    node_xsd = childNode
                    break
        if ( node_xsd == None):
            raise ValueError("Undefined type in schema: " + stType)

        #return f_xml_GetChildElement(node_xsd)[0] #Return complexType.sequence
        return node_xsd


    def Serilize(self,obj):
        '''
        Serilize object to a minidom xml document
        Objects' type must be defined in schema
        '''
        if ( self.schema == None ):
            raise ValueError("XMLSerilizer must be assigned with XML schema before use")

#Start parse object one by one for all element in schema
#Create new xml
        doc = minidom.Document()
#Find first element in xsd and treat it as the object 
        element1 = [x for x in f_xml_GetChildElement(self.schema) if x.nodeName=="xs:element"][0]
        doc.appendChild(self.__SerilizeNode__(obj,element1,doc))
        return doc
        
    def SerilizeToFile(self,obj,stFileName):
        '''
        Serilize object to a xml file
        Reference 'Serilize'.
        '''
        doc = self.Serilize(obj)
        f = open(stFileName,'w')
        doc.writexml(f,indent="",addindent="    ",newl="\n",encoding="UTF-8")
        f.close()
#reformat
        f = open(stFileName,"r")
        st = f.read()
        f.close()
        st = f_xml_TrimAllText(st)
        f = open(stFileName,'w')
        f.write(st)
        f.close()

        return doc
    
    def __append_child__(self,node,child):
        '''
        Wrapper of node.appendChild function'
        Additional feature: if the child is None, then do not do anything
        '''
        if (child is not None):
            node.appendChild(child)
        return

    def __SerilizeNode__(self,obj,node_xsd,doc):
        '''
        Serilize object with given node in xsd, and append it into xml tree
        Iterative
        '''
        #print("Seri: ",obj,node_parent.toxml(),node_xsd.toxml())
#Null object
        if ( obj == None ):
            stMin = node_xsd.getAttribute("minOccurs")
            if ( stMin == ""):
                stMin = "1"
            if ( int(stMin) > 0 ):
                raise ValueError("Null object which must be non-Null in schema: %s" % node_xsd.getAttribute("name"))
            #else: #Return empty text node, which is available to append
            #    return doc.createTextNode("")
            else: #Return none, which means do not append this
                return None

#       if ( node_xsd.nodeName == "xs:extension"):
#           stType = node_xsd.getAttribute("base")
#           stNodeName = "__base__"
#       else:
#           stType =  node_xsd.getAttribute("type")
#           stNodeName = node_xsd.getAttribute("name")
        stType =  node_xsd.getAttribute("type")
        stNodeName = node_xsd.getAttribute("name")

        #New element
        value = None
        if ( f_xml_IsBaseType(stType)):
            value = str(obj)
            #Deal with element / attribute
            if ( node_xsd.nodeName == "xs:attribute"):
                node_now = doc.createAttribute(stNodeName)
                node_now.nodeValue = value
            elif ( node_xsd.nodeName == "xs:element"):
                node_now = doc.createElement(stNodeName)
                node_text = doc.createTextNode(value)
                node_now.appendChild(node_text)
        else:
            node_now = doc.createElement(stNodeName)
#Complex type, iterate all of its elements
            #Look for type definition for element
            node_xsd_new = self.__FindTypeSchema__(stType)
#Array type, enumerate all of its content in order
            if ( self.__IsArray__(node_xsd_new)):
#For array, use element one by one
                nCount = 0
                for childNode in node_xsd_new.getElementsByTagName("xs:element"):
                    nStart = nCount
                    nEnd = nCount+1
#As we don't know the name in list, so we must present a fixed number of elements. Use maxOccurs.
                    stMax =childNode.getAttribute("maxOccurs").encode("ascii")  
                    if ( stMax == "unbounded" ):
                        nEnd = len(obj)
                    else:
                        nEnd = nStart + int(stMax)
                    #print("Iter :",nEnd)
                    for i in range(nStart,nEnd):
                        #print(i)
#                       node_now.appendChild(self.__SerilizeNode__(obj[i],childNode,doc))
                        self.__append_child__(node_now,self.__SerilizeNode__(obj[i],childNode,doc))
#Update position
                    nCount = nEnd
            else:
#Deal with base in extension first
                node_xsd_ext_type = node_xsd_new.getElementsByTagName("xs:extension")
                if ( len(node_xsd_ext_type) == 1):
                    node_xsd_ext_type = node_xsd_ext_type[0]
                    #Reset definition to base, deal with it and set back for later 
                    node_xsd.setAttribute("type", node_xsd_ext_type.getAttribute("base"))
                    if ( hasattr(obj,"__base__")):
                        node_now = self.__SerilizeNode__(getattr(obj,"__base__"),node_xsd,doc) 
                    else:
                        node_now = self.__SerilizeNode__(obj,node_xsd,doc)
                    node_xsd.setAttribute("type", node_xsd_new.getAttribute("name"))

#Parse object
#Attribute first, element second
                listChild = node_xsd_new.getElementsByTagName("xs:attribute") +node_xsd_new.getElementsByTagName("xs:element") 

                if ( isinstance(obj,list)):
                    #Its children is stored in list
                    nCount = 0
                    for childNode in listChild:
                        nStart = nCount
                        nEnd = nCount+1
#One attribute always has one value
                        if ( childNode.nodeName == "xs:attribute"):
                            stMax = 1
                        else:
#As we don't know the name in list, so we must present a fixed number of elements. Use maxOccurs.
                            stMax =childNode.getAttribute("maxOccurs").encode("ascii")  
                        if ( stMax == "unbounded" ):
                            nEnd = len(obj)
                        else:
                            nEnd = nStart + int(stMax)
                        #print("Iter :",nEnd)
                        for i in range(nStart,nEnd):
                            #print(i)
#Array must be element so we do not need to deal with attribute
                            self.__append_child__(node_now,self.__SerilizeNode__(obj[i],childNode,doc)) 
                           #node_now.appendChild(self.__SerilizeNode__(obj[i],childNode,doc))
                else:
                    for childNode in listChild:
                        name1 = childNode.getAttribute("name") 
                        childObj = None
                        if (hasattr(obj,name1)):
                            childObj = getattr(obj,name1)
                        #Detect list or single element
                        if ( self.__IsArrayElement__(childNode)):
#Array must be element so we do not need to deal with attribute
                            for childObj2 in childObj:
#                               node_now.appendChild(self.__SerilizeNode__(childObj2,childNode,doc))
                                self.__append_child__(node_now,self.__SerilizeNode__(childObj2,childNode,doc))
                        else:
#Add by element/attribute
                            if ( childNode.nodeName == "xs:element"):
#                                node_now.appendChild(self.__SerilizeNode__(childObj,childNode,doc))
                                self.__append_child__(node_now,self.__SerilizeNode__(childObj,childNode,doc))
                            elif ( childNode.nodeName == "xs:attribute"):
                                node_attr = self.__SerilizeNode__(childObj,childNode,doc) 
                                node_now.setAttributeNode(node_attr)
                            else:
                                raise ValueError("Unknown nodeName %s" % childNode.nodeName)


            #Change normal xml element to text node if type is tablerow
            if ( self.__IsTable__(node_xsd_new)):
                stNew = ""
#Clear tag in tablerow
                if ( "TableRow" in node_xsd_new.getAttribute("name")):
                    for i,childNode in enumerate(node_now.childNodes):
                        #print(childNode.toxml())
                        if ( i != 0):
                            stNew += " "
#join all text node data in xml tree
                        stNew += childNode.childNodes[0].data
#Clear tag in table
                else:
                    for i,childNode in enumerate(node_now.childNodes):
                        if ( i != 0):
                            stNew += "\n"
                        stNew += childNode.childNodes[0].data
                    #Append new line before and after table (use empty text node instead)
#                   stNew = "".join(["\n",stNew,"\n"])

                if ( node_xsd.nodeName == "xs:attribute"):
                    node_now = doc.createAttribute(stNodeName)
                    node_now.nodeValue = stNew 
                elif ( node_xsd.nodeName == "xs:element"):
                    node_now = doc.createElement(stNodeName)
                    node_text = doc.createTextNode(stNew)
                    node_now.appendChild(node_text)
                    if ( "TableRow" in node_xsd_new.getAttribute("name")):
                        node_now.appendChild(node_text)
                    else:
                        node_now.appendChild(doc.createTextNode(""))
                        node_now.appendChild(node_text)
                        node_now.appendChild(doc.createTextNode(""))

        #print("Ser add:",node_now.toxml())
        return node_now


    def Deserilize(self,stFileName):
        '''
        Deserilize a xml file and return a new object
        The object must have a 0-parameter construction method
        '''
        doc = minidom.parse(stFileName)
#Root element is the basic type
        node = f_xml_GetChildElement(doc)[0]
#Create root object
#        obj = locals()[node.nodeName]()
#Find root type
        node_xsd = filter(lambda x: x.nodeName=="xs:element",f_xml_GetChildElement(self.schema))[0]
        #__FindTypeSchema__(node.nodeName) 
        obj = self.__DeserilizeNode__(node,doc,node_xsd)
        return obj

    def __DeserilizeNode__(self,node_xml,doc,node_xsd):
        '''
        Deserilize a xml node as a object's element
        Iterative
        '''
        #print("Deser pos:",node_xml.toxml(),node_xsd.toxml())
        #print("Deser pos:",node_xsd.getAttribute("name"))
        #Create main object
        obj = None
        stType = node_xsd.getAttribute("type").encode("ascii")
        if ( f_xml_IsBaseType(stType)):
            obj = f_xml_ReadBaseTypeNode(node_xml,stType)
            if ( obj == None): #Not found, try to use default
                if ( node_xsd.hasAttribute("default")):
                    obj = f_xml_ReadBaseTypeString(node_xsd.getAttribute("default"),stType)
                else:
                    raise ValueError("No value for node %s" % node_xml.toxml())
            #print("Create object: ",obj)
            return obj #Directly return if it is a simple type
        else:
            if ( self.__IsArray__(node_xsd)):
                obj = []
            else:
                if ( "TableRowOf" in stType):
                    obj = self.classDef[stType[10:]]()
                else:
                    obj = self.classDef[stType]()
            #print("Create object: ",obj)
        
#If it has child nodes then we proceed

#Lookup type definition
        node_type = self.__FindTypeSchema__(stType.encode("ascii"))
        node_mc_type = f_xml_GetChildElement(node_type)[0] #Get the main child, xs:sequence or xs:simpleContent. We assume the first one is the main node

#If it is Table/TableRow, we cannot use childNodes and corresponding nodes name, so we must dispatch XSD child and column indivdually
        if ( self.__IsTable__(node_xsd) ):
            if ( node_xsd.nodeName == "xs:element"):
                st_xml = node_xml.childNodes[0].data
            elif ( node_xsd.nodeName == "xs:attribute"):
                st_xml = node_xml.getAttribute(node_xsd.getAttribute("name"))

            obj = self.__DeserilizeTable__(st_xml,doc,node_xsd)
            return obj

#Look for all child in xsd and create all element ( set None if not found )
#Extension, Attribute and Element all considered, order as Ext/Attr/Ele

        listChild = []
        tag_element = 0
        tag_attr = 1
        tag_text = 2

        #Read attribute
#       for childNodeXSD in f_xml_GetChildElement(node_type):
#           if ( childNodeXSD.nodeName == "xs:attribute"):
#               listChild.append(childNodeXSD)
#       #Read element
#       if ( node_mc_type.nodeName == "xs:sequence"):
#           for childNodeXSD in f_xml_GetChildElement(node_mc_type):
#               if ( childNodeXSD.nodeName != "xs:element"):
#                   continue
#               listChild.append(childNodeXSD)

        listChild = node_type.getElementsByTagName("xs:attribute") +node_type.getElementsByTagName("xs:element") 

        if ( node_mc_type.nodeName == "xs:simpleContent"):
            node_ext_type = f_xml_GetChildElement(node_mc_type)[0]
            if ( node_ext_type.nodeName != "xs:extension"):
                raise ValueError("simpleContent must have one and only one element named extension")
            #Read base text
#            listChild.append(node_ext_type)

            #Invoke procedure to create base type of extension
            stExtType = node_ext_type.getAttribute("base")
            if ( f_xml_IsBaseType(stExtType)):
                obj.__base__ = f_xml_ReadBaseTypeNode(node_xml,stExtType)
            else:
            #Create a temporary element in xsd
                node_extbase = doc.createElement("xs:element")
                node_extbase.setAttribute("type",stExtType)
                obj2 = self.__DeserilizeNode__(node_xml,doc,node_extbase)
            #list and other basic types must be put into new obj
            #class directly expanded
                if ( isinstance(obj2,(list,str,int,float) )):
                    obj.__base__ = obj2
                else:
                    obj = obj2


#End parse xsd, start read xml
        listObj = [] #Store name and value
        tag_multi = 2
        tag_single = 1

        for childNodeXSD in listChild:
#Normal xml
#2 Condition : array or class
            #
#            if ( childNodeXSD.nodeName == "xs:extension"):
#                stChildType = childNodeXSD.getAttribute("base")
#                if ( f_xml_IsBaseType(childNodeXS
#                listObj.append(["__base__",f_xml_ReadBaseTypeString(node_xml.childNodes[0].data,childNodeXSD.getAttribute("base")),tag_single])
#                continue

            stNodeName = childNodeXSD.getAttribute("name")
            stChildType = childNodeXSD.getAttribute("type")
            #print("Type of child: " + stChildType)

            if ( childNodeXSD.nodeName == "xs:attribute"):
                st_xml = node_xml.getAttribute(stNodeName)
                if ( not f_xml_IsBaseType(stChildType)):
#Parse tablerow in attribute
                    if ( "TableRow" in stChildType ):
#Create an temporary element for table parsing
                        listObj.append([stNodeName,self.__DeserilizeTable__(st_xml,doc,childNodeXSD),tag_single])
                    else:
                        raise ValueError("Attribute must be a base type or table/tablerow!")
                else:
                    listObj.append([stNodeName,f_xml_ReadBaseTypeString(st_xml,stChildType),tag_single])
                continue
              
            
            listNodeXML = filter(lambda x:x.nodeName==stNodeName,f_xml_GetChildElement(node_xml))

            #Deal with cases when not found
            if ( len(listNodeXML) == 0):#Does not found,use None if no default
#If there is default value, create temporary node for processing
                if ( childNodeXSD.hasAttribute("default")):
                    node_default = doc.createElement(stNodeName)
                    node_text = doc.createTextNode(childNodeXSD.getAttribute("default"))
                    node_default.appendChild(node_text)
                    listNodeXML.append(node_default)
#Exception when it cannot be not none
                if ( int(childNodeXSD.getAttribute("minOccurs")) > 0 ):
                    raise ValueError("Cannot find element which must be non-Null in schema: " + stNodeName)

                if ( not isinstance(obj,list)):
                    setattr(obj,stNodeName,None)
                    #print("Set %s to None" % stNodeName)
                else:
#We do not deal with list here as elements in list can never be neglected
                    raise ValueError("Wrong XSD: Elements in list cannot be nulliable: %s" % stNodeName)

            #print("list node: ",[x.toxml() for x in listNodeXML])
            for childNodeXML in listNodeXML:
                obj2 = self.__DeserilizeNode__(childNodeXML,doc,childNodeXSD)
            #Deal with multiple one inlist as an attribute of class
                if ( childNodeXSD.getAttribute("maxOccurs") != "1" ):
                    listObj.append([stNodeName,obj2,tag_multi])
                else:
            #Deal with single attribute of class
                    listObj.append([stNodeName,obj2,tag_single])

#All child values are now stored, save them to parent object
        for kvp in listObj:
            stNodeName = kvp[0]
            obj2 = kvp[1]
            tag = kvp[2]
        #Deal with list parent
            if ( isinstance(obj,list)):
                obj.append(obj2)
            #Deal with class parent
            elif ( tag == tag_multi ):
                if( hasattr(obj,stNodeName)):
                    getattr(obj,stNodeName).append(obj2)
                else:
                    setattr(obj,stNodeName,[obj2]);
            else:
                setattr(obj,stNodeName,obj2)

        #print("Final return:",obj)
        return obj

    def __DeserilizeTable__(self,st_xml,doc,node_xsd):
        '''
        Deserilize a string to TableOf and TableRowOf types
        '''
        stType = node_xsd.getAttribute("type")
#Create object
        if ( self.__IsArray__(node_xsd)):
            obj = []
        else:
            if ( "TableRowOf" in stType):
                obj = self.classDef[stType[10:]]()
            else:
                raise ValueError("Unknown type %s" %stType)

#Lookup type definition
        node_type = self.__FindTypeSchema__(stType.encode("ascii"))
        node_mc_type = f_xml_GetChildElement(node_type)[0] #Get the main child, xs:sequence

#Table or TableRow condition
#Table  condition
        if ( "TableOf" in stType):
#Split by endofline
            listNodeSub = st_xml.strip().split("\n")
            rowNodes = f_xml_GetChildElement(node_mc_type)
#Count nodes by count(maxOccurs)
            nCount = 0
            for rowNode in rowNodes:
                nStart = nCount
                nEnd = nCount+1
#As we don't know the name in list, so we must present a fixed number of elements. Use maxOccurs.
                stMax =rowNode.getAttribute("maxOccurs").encode("ascii")  
                if ( stMax == "unbounded" ):
                    nEnd = len(listNodeSub)
                else:
                    nEnd = nStart + int(stMax)

                for i in range(nStart,nEnd):
#Create a temporary element for TableRow type to deal with
                    node_row = doc.createElement("Item")
                    node_row_text = doc.createTextNode(listNodeSub[i])
                    node_row.appendChild(node_row_text)
                    obj.append(self.__DeserilizeNode__(node_row,doc,rowNode))

#Table Row condition
#This is special as the input maybe a raw string (passed by table dispatcher), not xml node
        else:
            #print("Node_xml: " + node_xml)
            listNodeSub = st_xml.split()
            #print(listNodeSub)

            #Split by everything can be split
            columnNodes = f_xml_GetChildElement(node_mc_type)
            #Calculate column count
            listCol = []#Record column count
            listMulti = []#Record whether to use list
            nCol = 0
            for childNode in columnNodes:
                stMax =childNode.getAttribute("maxOccurs").encode("ascii")
                if ( stMax == "unbounded"):
                    if( len(listCol) != len(columnNodes)-1 ):
                         raise ValueError("only the last element in definition can be unbounded")
                    else:
                        nCol = sum(listCol)
                        if ( nCol > len(listNodeSub)):
                            raise ValueError("Too few columns to deserilize")
                        listCol.append(len(listNodeSub)-nCol)
                        listMulti.append(True)
                else:
                    listCol.append(int(stMax))
                    listMulti.append(listCol[-1] > 1)

            nCol = sum(listCol)
            if ( nCol != len(listNodeSub)):
                raise ValueError("Number of columns is different in object and definition")

            #Parse the string to object
            index = -1
            for j,nCol in enumerate(listCol):
                stNodeName = columnNodes[j].getAttribute("name")
                for i in range(0,nCol):
                    index += 1
                    #Create xml element to call __DeserilizeNode__
                    node_col = doc.createElement(stNodeName)
                    node_col_text = doc.createTextNode(listNodeSub[index])
                    node_col.appendChild(node_col_text)

                    #Parse
                    obj2 = self.__DeserilizeNode__(node_col,doc,columnNodes[j])
                    #Append
                    if ( isinstance(obj,list)):
                        obj.append(obj2)
                    #Deal with class parent
                    elif ( listMulti[j] ):
                        if( hasattr(obj,stNodeName)):
                            getattr(obj,stNodeName).append(obj2)
                        else:
                            setattr(obj,stNodeName,[obj2]);
                    else:
                        setattr(obj,stNodeName,obj2)

        #print("Final Table return:",obj)
        return obj

    @classmethod
    def load_xml(cls,vars2,cls2,filename):
        '''
        This function accept a class with variable filename_xsd,
        and return an object according to a XML file

        Note, because it is required all possible classes included present, so it is necessary to pass outside

        :param vars2: generally pass global() where it is used in the module
        '''
        if (hasattr(cls2,"filename_xsd")):
            filename_xsd = cls2.filename_xsd
            xp = cls(vars2,cls2.filename_xsd)
            obj = xp.Deserilize(filename)
            return obj
        else:
            raise ValueError("load_xml only works with classes with a class variable \"filename_xsd\"")

    @classmethod
    def save_xml(cls,cls2,obj,filename):
        '''
        Save parameters to a XML file

        :param cls2: the class of the object (We do not use type(obj) to avoid some problems because it is unreliable for old-style class)
        :param obj: the object
        '''
        if (hasattr(cls2,"filename_xsd")):
            filename_xsd = cls2.filename_xsd
            xp = cls(globals(),cls2.filename_xsd)
            xp.SerilizeToFile(obj,filename)
            return
        else:
            raise ValueError("load_xml only works with classes with a class variable \"filename_xsd\"")

