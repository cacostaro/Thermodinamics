import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from tabulate import tabulate
from scipy.interpolate import interp1d
from scipy.interpolate import UnivariateSpline
from scipy.optimize import curve_fit
#from scipy.signal import savitzky_golay


#############################              Script by Carlos Acosta              #############################



#####################################               Data               ######################################

R = 0
sustance = 'Isobutano'


w = 0.181
TcK = 408.1 
PcBar = 36.48
tfus = 113.8

def generate_Temperature(tcritic, n, tfus):
    
    sum = tfus
    delaTr = tcritic - tfus
    fact = delaTr / n

    temperatures = []

    temperatures_tYPE = [
                        113.8, 120, 140, 160, 180, 200, 220, 240, 260, 270, 280, 290, 300, 
                        310, 320, 330, 340, 350, 360, 370, 380, 390, 400, 408.1
                         ]
    '''
    while (sum < tcritic):

        temperatures.append(round(sum, 2))
        # Calcular el nuevo valor de la suma sumando el valor de "fact"
        sum += fact
        # Verificar si la suma supera 407
        if sum >= (tcritic):
            break  # Salir del bucle si la suma supera 407
    '''
    return temperatures_tYPE

def generate_Pressure(pcritic, n):
    
    sum = 1
    
    fact = pcritic / n

    pressures = []

    

    pressures_TYPE = [0.00000019, 0.00000093, 0.000048, 0.00082, 0.007, 0.0369, 0.1374, 0.3989, 0.9600, 1.4081,
                     2.0020, 2.7686, 3.7365, 4.9340, 6.3920, 8.1400, 10.2100, 12.6400, 15.4600, 18.7200, 22.4800, 
                     26.8200, 31.8600, 36.4976
                     ]
    



    while (sum < pcritic):

        pressures.append(round(sum, 2))
        # Calcular el nuevo valor de la suma sumando el valor de "fact"
        sum += fact
        # Verificar si la suma supera 407
        if sum >= (pcritic-0.5):
            break  # Salir del bucle si la suma supera 407

    return pressures_TYPE

###################################             NO  Modify               ####################################
print("__________________________________________________________________\n")
print("======= SOLUCION ECUACIONES CUBICAS DE ESTADO  PARA", sustance.upper() ,"=======\n")
num_Data = int(input("Ingrese cauntos valores para temperatura y presion desea: "))
print(generate_Temperature(TcK, num_Data, tfus))
print(generate_Pressure(PcBar, num_Data))


# Obtener los datos de temperatura y presión
temperatures = generate_Temperature(TcK, num_Data, tfus)
pressures = generate_Pressure(PcBar, num_Data)



# Asegurarse de que ambas listas tengan la misma longitud


# Combinar las listas de temperaturas y presiones
combined_data = [[temp, pres] for temp, pres in zip(temperatures, pressures)]

# Imprimir la tabla combinada
print("Tabla combinada de Temperatura y Presion:")
print(tabulate(combined_data, headers=["Temperatura (K)", "Presion (Bar)"], tablefmt="fancy_grid", showindex=False))
#print(tabulate(map(lambda x: [x[0], x[1]], combined_data), headers=["Temperatura (K)", "Presion"], tablefmt="fancy_grid", showindex=False))



# Crear un DataFrame de pandas con los datos de la tabla
#df = pd.DataFrame(combined_data, columns=["Temperatura (K)", "Presion"])

# Guardar el DataFrame como un archivo Excel
#nombre_archivo = "datos_temperatura_presion.xlsx"
#df.to_excel(nombre_archivo, index=False)

#print(f"Los datos se han exportado correctamente a '{nombre_archivo}'")



######################################              SRK               #######################################

def calculate_Tr_values(temperatures, tcritic):
    Tr_values = []
    for temp in temperatures:
        Tr = temp / tcritic
        Tr_values.append(round(Tr, 7))
    return Tr_values

def calculate_Pr_values(pressures, pcritic):
    Pr_values = []
    for pres in pressures:
        Pr = pres / pcritic
        Pr_values.append(round(Pr, 7))
    return Pr_values

def calculate_Alfa_values(Tr_values, w):
    Alfa_values = []
    for Tr in Tr_values:
        Alfa = (1+(0.48+(1.475*w)-0.176*w**2)*(1-Tr**0.5))**2
        Alfa_values.append(round(Alfa, 7))
    return Alfa_values

def calculate_A_values(Alfa_values, Pr_values, Tr_values):
    A_values = []
    for Alfa, Pr, Tr in zip(Alfa_values, Pr_values, Tr_values):
        A = 0.42748 * (Alfa * Pr / Tr**2)
        A_values.append(round(A, 8))
    return A_values

def calculate_B_values(Pr_values, Tr_values):
    B_values = []
    for Pr, Tr in zip(Pr_values, Tr_values):
        B = 0.08664*(Pr/Tr)
        B_values.append(round(B, 8))

    return B_values



def generate_table(temperatures, pressures, Tr_values, Pr_values, Alfa_values, A_values, B_values):
    # Combinar los valores en una lista de listas
    combined_data = []
    for temp, pres, Tr, Pr, Alfa, A,  B in zip(temperatures, pressures, Tr_values, Pr_values, Alfa_values, A_values, B_values):
        combined_data.append([temp, pres, Tr, Pr, Alfa, A, B])

    # Imprimir la tabla
    print("Tabla combinada de Temperatura y Presion:")
    print(tabulate(combined_data, headers=["Temperatura (K)", "Presion (Bar)", "Tr", "Pr", "Alpha", "A", "B"], tablefmt="fancy_grid", showindex=False))

Tr_values = calculate_Tr_values(temperatures, TcK)
Pr_values = calculate_Pr_values(pressures, PcBar)
Alfa_values = calculate_Alfa_values(Tr_values, w)
A_values = calculate_A_values(Alfa_values, Pr_values, Tr_values)
B_values = calculate_B_values(Pr_values, Tr_values)


# Generar la tabla con los valores adicionales
generate_table(temperatures, pressures, Tr_values, Pr_values, Alfa_values, A_values, B_values)




##################################            Newton Raphson              ##################################

def newton_raphson_met1(A, B, tolerancia=1e-8, max_iter=1000):
    # Funcion F(z) y su derivada F'(z)
    def F(Z, A, B):
        return Z ** 3 - Z ** 2 + Z * (A - B - B ** 2) - A * B

    def F_prime(Z, A, B):
        return 3 * Z ** 2 - 2 * Z + (A - B - 2 * B ** 2)

    # Estimacion inicial para Z
    Z = 0
    F_values = []
    F_prime_values = []


    # Proceso iterativo del metodo de Newton-Raphson
    for _ in range(max_iter):
        F_Z = F(Z, A, B)
        F_prime_Z = F_prime(Z, A, B)

        F_values.append(round(F_Z, 8))
        F_prime_values.append(round(F_prime_Z, 8))

        if abs(F_prime_Z) < tolerancia:
            break

        Z -= F_Z / F_prime_Z

    return round(Z, 8)#, F_values, F_prime_values


def newton_raphson_met2(A, B, tolerancia=1e-8, max_iter=1000):
    
    def F(Z, A, B):
        return Z ** 3 - Z ** 2 + Z * (A - B - B ** 2) - A * B

    def F_prime(Z, A, B):
        return 3 * Z ** 2 - 2 * Z + (A - B - 2 * B ** 2)

    
    Z = 1
    F_values = []
    F_prime_values = []


    
    for _ in range(max_iter):
        F_Z = F(Z, A, B)
        F_prime_Z = F_prime(Z, A, B)

        F_values.append(round(F_Z, 8))
        F_prime_values.append(round(F_prime_Z, 8))


        if abs(F_prime_Z) < tolerancia:
            break

        Z -= F_Z / F_prime_Z

    return round(Z, 8)#, F_values, F_prime_values





def generate_table_with_Z(temperatures, pressures, Tr_values, Pr_values, Alpha_values, A_values, B_values, Z_values):
    # Combinar los valores en una lista de listas
    combined_data = []
    for temp, pres, Tr, Pr, Alpha, A, B, Z in zip(temperatures, pressures, Tr_values, Pr_values, Alpha_values, A_values, B_values, Z_values):
        combined_data.append([temp, pres, Tr, Pr, Alpha, A, B, Z])

    # Imprimir la tabla
    headers = ["Temperatura (K)", "Presion (Bar)", "Tr", "Pr", "Alpha", "A", "B", "Z", "Volumen específico"]
    v_values = calculate_specific_Volume(Z_values, temperatures, temperatures)
    for data, v in zip(combined_data, v_values):
        data.append(v)
    print(tabulate(combined_data, headers=headers, tablefmt="fancy_grid", showindex=False))








def calculate_specific_Volume(Z_values, temperatures, pressures):
    # Constante de los gases
    R_gas = 0.08314472  # Bar L / (K mol)
    
    # Calcular el volumen específico para cada valor de Z
    v_values = []
    for Z, T, P in zip(Z_values, temperatures, pressures):
        v = (Z * R_gas * T) / P
        v_values.append(v)
    return v_values



#Z_values = [newton_raphson_method(A_values, B_values) for _ in range(num_Data)]   
#Z_values = [newton_raphson_method(A, B) for A, B in zip(A_values, B_values)]

Z1_values = [newton_raphson_met1(A, B) for A, B in zip(A_values, B_values)]
Z2_values = [newton_raphson_met2(A, B) for A, B in zip(A_values, B_values)]


specific_Volumez1 = calculate_specific_Volume(Z1_values, temperatures, pressures)[1:]   # OJOOOOO
specific_Volumez2 = calculate_specific_Volume(Z2_values, temperatures, pressures)[1:]



# Generar la tabla con los valores de Z
generate_table_with_Z(temperatures, pressures, Tr_values, Pr_values, Alfa_values, A_values, B_values, Z1_values)
generate_table_with_Z(temperatures, pressures, Tr_values, Pr_values, Alfa_values, A_values, B_values, Z2_values)



# Crear un DataFrame de pandas con los datos de la tabla
df1 = pd.DataFrame({
    "Temperatura (K)": temperatures,
    "Presion(Bar)": pressures,
    "Tr": Tr_values,
    "Pr": Pr_values,
    "Alfa": Alfa_values,
    "A": A_values,
    "B": B_values,
    "Z": Z1_values
})

df2 = pd.DataFrame({
    "Temperatura (K)": temperatures,
    "Presion(Bar)": pressures,
    "Tr": Tr_values,
    "Pr": Pr_values,
    "Alfa": Alfa_values,
    "A": A_values,
    "B": B_values,
    "Z": Z2_values
})

'''
# Guardar el DataFrame como un archivo Excel
name_Data = "datos_temperatura_presionTEST003.xlsx"
nombre_hoja1 = "Hoja1"  # Nombre de la primera hoja
nombre_hoja2 = "Hoja2"  # Nombre de la segunda hoja
with pd.ExcelWriter(name_Data) as writer:
    df1.to_excel(writer, sheet_name=nombre_hoja1, index=False)
    df2.to_excel(writer, sheet_name=nombre_hoja2, index=False)

print(f"Los datos se han exportado correctamente a '{name_Data}'")
'''


#\ Grafico

plt.plot(calculate_specific_Volume(Z1_values, temperatures, pressures), pressures, marker='o', linestyle='-')
plt.xlabel('Volumen Específico')
plt.ylabel('Temperatura (K)')
plt.title('Gráfico de Líneas: Volumen Específico vs. Presion')
plt.grid(True)
plt.show()

plt.plot(calculate_specific_Volume(Z2_values, temperatures, temperatures), pressures, marker='o', linestyle='-')
plt.xlabel('Volumen Específico')
plt.ylabel('Temperatura (K)')
plt.title('Gráfico de Líneas: Volumen Específico vs. Presion')
plt.grid(True)
plt.show()



# Gráfico de Z1_values
plt.plot(calculate_specific_Volume(Z1_values, temperatures, temperatures), pressures, marker='o', linestyle='-', label='Z1')

 #Gráfico de Z2_values
plt.plot(calculate_specific_Volume(Z2_values, temperatures, temperatures), pressures, marker='o', linestyle='-', label='Z2')

'''
#Isotermas
'''
plt.xlabel('Volumen Específico')
plt.ylabel('Presion (Bar)')
plt.title('Gráfico de Líneas: Volumen Específico vs. Presion')
plt.grid(True)
plt.legend()  # Mostrar leyenda con etiquetas Z1 y Z2
plt.show()




volumen_Z1 = calculate_specific_Volume(Z1_values, temperatures, temperatures)
volumen_Z2 = calculate_specific_Volume(Z2_values, temperatures, temperatures)


# Crea un DataFrame de pandas con los valores de volumen específico y presión
data = {'Presión': pressures,
        'Volumen Específico Z1': volumen_Z1,
        'Volumen Específico Z2': volumen_Z2,
        }

dfz = pd.DataFrame(data)

# Imprime el DataFrame
print("Tabla de valores:")
print(dfz)



