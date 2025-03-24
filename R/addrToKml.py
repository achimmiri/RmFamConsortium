import urllib
import urllib.request
import urllib.parse
import xml.dom.minidom
import json 
import pandas as pd

def geocode(address, sensor=False):
 # This function queries the Google Maps API geocoder with an
 # address. It gets back a csv file, which it then parses and
 # returns a string with the longitude and latitude of the address.

 # This isn't an actual maps key, you'll have to get one yourself.
 # Sign up for one here: https://code.google.com/apis/console/
  mapsKey = 'xxxx-xxxx-xxxx-xxxx'
  mapsUrl = 'https://maps.googleapis.com/maps/api/geocode/json?address='
     
 # This joins the parts of the URL together into one string.
  url = ''.join([mapsUrl,urllib.parse.quote(address),'&sensor=',str(sensor).lower(), '&key=',mapsKey])
  jsonOutput = str(urllib.request.urlopen(url).read()) # get the response
  # fix the output so that the json.loads function will handle it correctly
  jsonOutput=jsonOutput.replace ("\\n", " ")
  jsonOutput = jsonOutput.replace("\\'", " ")
  jsonOutput = jsonOutput.replace("\\", "\\\\")
  print(jsonOutput)
  result = json.loads(jsonOutput[2:-1]) # converts jsonOutput into a dictionary

  # check status is ok i.e. we have results (don't want to get exceptions)
  if result['status'] != "OK": 
    return ""
  coordinates=result['results'][0]['geometry']['location'] # extract the geometry 
  return str(coordinates['lng'])+','+str(coordinates['lat'])

def createKML(places, fileName):
 # This function creates an XML document and adds the necessary
 # KML elements.

  kmlDoc = xml.dom.minidom.Document()
  
  kmlElement = kmlDoc.createElementNS('http://earth.google.com/kml/2.2','kml')

  kmlElement = kmlDoc.appendChild(kmlElement)

  documentElement = kmlDoc.createElement('Document')
  documentElement = kmlElement.appendChild(documentElement)

  for address in places:
   placemarkElement = kmlDoc.createElement('Placemark')

   descriptionElement = kmlDoc.createElement('description')
   descriptionText = kmlDoc.createTextNode(address)
   descriptionElement.appendChild(descriptionText)
   placemarkElement.appendChild(descriptionElement)
   pointElement = kmlDoc.createElement('Point')
   placemarkElement.appendChild(pointElement)
   coorElement = kmlDoc.createElement('coordinates')

   # This geocodes the address and adds it to a  element.
   coordinates = geocode(address)
   coorElement.appendChild(kmlDoc.createTextNode(coordinates))
   pointElement.appendChild(coorElement)

   documentElement.appendChild(placemarkElement)

  # This writes the KML Document to a file.
  kmlFile = open(fileName, 'w')
  kmlFile.write(kmlDoc.toprettyxml(' '))  
  kmlFile.close()

if __name__ == '__main__':
  wb = pd.read_excel("Contributor-Affliciations.xlsx")
  # addr=[]
  # for x in wb['Affliciation ']:
  #  if not isinstance(x,(float)):
  #   addr.append(x.replace("\n"," "))

  addr = [x.replace("\n","") for x in wb['Affliciation '] if not isinstance(x, float)]

  createKML(addr, 'google.kml')
