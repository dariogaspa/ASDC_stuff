from xml.dom import minidom
xmldoc = minidom.parse('0902_fit1.xml')
itemlist = xmldoc.getElementsByTagName('source')
reg=open("file.reg",'w')

print len(itemlist)
print itemlist[0].attributes['name'].value
for s in itemlist :
    
    namesource= s.attributes['name'].value
#for node in xmldoc.getElementsByTagName('source'):
   #print node.toxml()
 #   print node.childNodes
    
    #print s.childNodes[3].childNodes[1].getAttribute("name")
    ravalue=s.childNodes[3].childNodes[1].getAttribute("value")
    #print s.childNodes[3].childNodes[3].getAttribute("name")
    decvalue=s.childNodes[3].childNodes[3].getAttribute("value")
    reg.write("fk5;point(%s,%s)# point=cross  text={%s} \n" % (ravalue, decvalue, namesource))
