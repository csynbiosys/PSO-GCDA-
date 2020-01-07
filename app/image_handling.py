import os 
from PIL import Image, ImageFile
from io import BytesIO
import base64
import sys
from scipy.interpolate import UnivariateSpline
from numpy import linspace, array
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pickle import dump
import traces

def decode_io(circuit_name, io):
    pngtxt = open('app/static/' + io + '/' + circuit_name, 'rb').read()
    print(pngtxt)
    im = Image.open(BytesIO(base64.b64decode(pngtxt)))
    im.save('app/static/' + io + '/' + circuit_name, 'PNG')
    os.remove('app/static/' + io + '/' + circuit_name + '.txt') 
    
def rgb2hex(r, g, b):
    return '#{:02x}{:02x}{:02x}'.format(r, g, b)

def hex_to_int(s):
    return [int(s[i:i+2], 16) for i in range(1, 7, 2)]

def hexDifference(hex1, hex2):
    return sum(abs(i - j) for i, j in zip(*map(hex_to_int, (hex1, hex2))))

def color_name(hexCol, colorDict):
    return min([(n, hexDifference(hexCol, c)) for n, c in colorDict.items()], key=lambda t: t[1])[0]
    
def extract_time_series(circuit_name, nodes, folder):
    image = Image.open('app/static/' + str(folder) +'/' + str(circuit_name))
    
    if image.mode in ('RGBA', 'LA') or (image.mode == 'P' and 'transparency' in image.info):
        pixels = image.convert('RGBA').load()
    
    xyHex = []
    for x in range(750):
        for y in range(450):
            r, g, b, a = pixels[x, y]
            hexCol = rgb2hex(r, g, b)
            if hexCol != '#000000': #filter background
                xyHex.append((x, y, hexCol))
                
    colorDict = {'purple': '#ae00ff', 
                 'dark_green': '#265300', 
                 'light_green': '#75ff00', 
                 'yellow': '#ffcf33', 
                 'brown': '#986928', 
                 'light_blue': "#00f8ff", 
                 'dark_blue': '#110e7a', 
                 'red': '#ff0900', 
                 'black': '#000000', 
                 'orange': '#ffab00', 
                 'pink': '#ff008e', 
                 'grey': '#898989'}
    
    sorted_color_Dict = {'purple': [], 
                         'dark_green': [], 
                         'light_green': [], 
                         'yellow': [], 
                         'brown': [], 
                         'light_blue': [], 
                         'dark_blue': [], 
                         'red': [], 
                         'black': [], 
                         'orange': [], 
                         'pink': [], 
                         'grey': []}
    
    for col in xyHex:
        sorted_color_Dict[color_name(col[2], colorDict)].append((col[0], col[1]))
    
    real_time_series = sorted([(k, len(sorted_color_Dict[k])) for k in sorted_color_Dict.keys()], key=lambda x: x[1])[-int(nodes):]
    
    real_time_series_keys = [t[0] for t in real_time_series]
    
    image_time_series = [sorted_color_Dict[k] for k in real_time_series_keys]
    
    return image_time_series

def rescale(image_time_series, max_expression, time_span):
    rescale_x = lambda x: (x/749.) * float(time_span)
    rescale_y = lambda y: ((abs(y - 449.))/449) * float(max_expression)
    return [(rescale_x(pix[0]), rescale_y(pix[1])) for pix in image_time_series]

def smoothing(time_series, time_span):
    filtered_time_series = list(dict(time_series).items())
    uspl = UnivariateSpline([x[0] for x in filtered_time_series], [y[1] for y in filtered_time_series])
    time = linspace(0, int(time_span), int(time_span))
    return uspl(time)

def plot(smoothed_time_series, time_span, count, circuit_name, max_expression, folder):
    f = open('app/static/' + folder + '/extracted/' + circuit_name + str(count) + '.pickle', 'wb')
    dump(smoothed_time_series, f)
    f.close()
    plt.figure()
    plt.axis([0, int(time_span), 0, int(max_expression)])
    plt.plot(list(range(int(time_span))), smoothed_time_series)
    plt.xlabel('time(min)')
    plt.ylabel('protein (arbitrary units)')
    plt.savefig('app/static/' + folder +'/extracted/' + circuit_name + str(count) + '.png')
    
def image_handling(circuit_name, max_expression, time_span, nodes, folder):
    ets = extract_time_series(circuit_name, nodes, folder)
    rescaled = [rescale(ts, max_expression, time_span) for ts in ets]
    smoothed = [smoothing(ts, time_span) for ts in rescaled]
    print(rescaled, file=sys.stderr)
    [plot(s_ts[1], time_span, s_ts[0], circuit_name, max_expression, folder) for s_ts in enumerate(smoothed)]
