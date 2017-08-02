import decimal as dc
import fractions as frac
import numpy as np
#%%
def fracspace(start, stop, size, incluir_start=True):
    """Puede no incluir el valor stop exacto, sino uno ligeramente menor."""
    step = (stop-start)/(size-1) if incluir_start==True else (stop-start)/(size)
    with dc.localcontext() as c:
        c.prec     = 5
        c.rounding = dc.ROUND_FLOOR
        step_dec = dc.Decimal(step)
    rango = range(size) if incluir_start==True else range(1, size+1)
    fracs = [frac.Fraction.from_decimal(start + x * step_dec) for x in rango]
    nums = np.array([x.numerator for x in fracs])
    denoms = np.array([x.denominator for x in fracs])
    return(nums, denoms)
#%%
# fracspace(1,10,18, incluir_start=False)
# fracspace(0,20,21)

# En caso de que prefiera limitar el tamaño de numerador y denominador, debería
# usar el método limit_denominator()
# Lo hago en la función fracspace2():

def fracspace2(start, stop, size, incluir_start=True, desarmados=False):
    """Puede no incluir el valor stop exacto, sino uno ligeramente menor. Limita
    el tamaño del denominador a 10 veces el valor del argumento size."""
    step = (stop-start)/(size-1) if incluir_start==True else (stop-start)/(size)
    with dc.localcontext() as c:
        c.prec     = 5
        c.rounding = dc.ROUND_FLOOR
        step_dec = dc.Decimal(step)
    rango = range(size) if incluir_start==True else range(1, size+1)
    fracs = [frac.Fraction.from_decimal(start + x * step_dec).limit_denominator(10*size) for x in rango]
    nums = ([x.numerator for x in fracs])
    denoms = ([x.denominator for x in fracs])
    return((nums, denoms) if desarmados==True else fracs)
#%%

#Ejemplos de prueba

# fracspace2(0,10,13, desarmados=False)
# fracspace2(1,10,100000)[:10]
# fracspace(1,10,100000)[:10]

#%%
# Supongamos que queremos generar un gráfico de P(Z=z) para valores de z entre
# 0 y 35, donde Z=X/Y y X,Y son poissonianas, entonces necesitamos exportar
# las siguientes listas para importarlas en Mathematica:

nums, denoms = fracspace2(0, 35, 500, incluir_start=False, desarmados=True)

np.savetxt('nums_fracspace2(0,35,500,incluir_start_False).csv', nums, delimiter=',')
np.savetxt('denoms_fracspace2(0,35,500,incluir_start_False).csv', denoms, delimiter=',')
