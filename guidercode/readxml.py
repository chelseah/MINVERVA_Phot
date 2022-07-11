import xml.etree.ElementTree as ET

def xmltodic(xmldata):
    data = ET.parse(xmldata)
    root = data.getroot()

    outdic = {}
    for i in range(len(root)):
        name1 = root[i].tag
        for sts in root[i]:
            name2 = sts.tag
            value = sts.text

            outdic['%s_%s' % (name1,name2)] = value

    return outdic

if __name__ == "__main__":
    import sys
    dic = xmltodic(sys.argv[1])

    print(dic['mount_azm_radian'])

    
