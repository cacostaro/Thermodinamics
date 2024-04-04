import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from tabulate import tabulate


#############################              Script by Carlos Acosta              #############################



#####################################               Data               ######################################

R = 0
sustance = 'Isobutano'


w = 0.181
TcK = 408.1 + 10
PcBar = 36.48
tfus = 142

def generate_Temperature(tcritic, n, tfus):
    
    sum = tfus
    delaTr = tcritic - tfus
    fact = delaTr / n

    temperatures = []

    while (sum < tcritic):

        temperatures.append(round(sum, 2))
        # Calcular el nuevo valor de la suma sumando el valor de "fact"
        sum += fact
        # Verificar si la suma supera 407
        if sum >= (tcritic-2):
            break  # Salir del bucle si la suma supera 407

    return temperatures

def generate_Pressure(pcritic, n):
    
    sum = 1
    
    fact = pcritic / n

    pressures = []

    while (sum < pcritic):

        pressures.append(round(sum, 2))
        # Calcular el nuevo valor de la suma sumando el valor de "fact"
        sum += fact
        # Verificar si la suma supera 407
        if sum >= (pcritic-0.5):
            break  # Salir del bucle si la suma supera 407

    return pressures

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

def newton_raphson_method(A_values, B_values):
    # Función F(z) y su derivada F'(z)
    def F(Z, A, B):
        return (Z**3) - (Z**2) + (Z * (A - B - B**2)) - (A * B)

    def F_prime(Z, A, B):
        return (3 * Z**2) - (2*Z) + (A - B - B**2)

    # Proceso iterativo del metodo de Newton-Raphson
    Z = 1 # Idealidad
    while True:
        # Evaluar F(Z) y su derivada en Z
        F_Z = F(Z, A_values, B_values)
        F_prime_Z = F_prime(Z, A_values, B_values)

        # Verificar si la derivada se esta acercando a cero
        if abs(F_prime_Z) < 1e-8:
            # La derivada se está acercando a cero, detener la funcion
            break

        # Actualizar Z utilizando la formula del método de Newton-Raphson
        Z_prime = Z - (F_Z / F_prime_Z)

        # Verificar convergencia (si Z no cambia o la diferencia es menor que 1e-8)
        if round(Z_prime, 8) == round(Z, 8):
            return round(Z_prime, 8)  # Devolver el valor de Z con hasta 8 

        # Actualizar Z para la proxima iteracion...
        Z = Z_prime

    return round(Z, 8)  # Devolver el valor de Z con hasta 8






def generate_table_with_Z(temperatures, pressures, Tr_values, Pr_values, Alpha_values, A_values, B_values, Z_values):
    # Combinar los valores en una lista de listas
    combined_data = []
    for temp, pres, Tr, Pr, Alpha, A, B, Z in zip(temperatures, pressures, Tr_values, Pr_values, Alpha_values, A_values, B_values, Z_values):
        combined_data.append([temp, pres, Tr, Pr, Alpha, A, B, Z])

    # Imprimir la tabla
    headers = ["Temperatura (K)", "Presion (Bar)", "Tr", "Pr", "Alpha", "A", "B", "Z", "Volumen específico"]
    v_values = calculate_specific_Volume(Z_values, temperatures, pressures)
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

Z_values = [newton_raphson_method(A, B) for A, B in zip(A_values, B_values)]

specific_Volume = calculate_specific_Volume(Z_values, temperatures, pressures)[1:]

# Generar la tabla con los valores de Z
generate_table_with_Z(temperatures, pressures, Tr_values, Pr_values, Alfa_values, A_values, B_values, Z_values)



# Grafico

plt.plot(calculate_specific_Volume(Z_values, temperatures, temperatures)[1:], pressures[1:], marker='o', linestyle='-')
plt.xlabel('Volumen Específico')
plt.ylabel('Temperatura (K)')
plt.title('Gráfico de Líneas: Volumen Específico vs. Temperatura')
plt.grid(True)
plt.show()




