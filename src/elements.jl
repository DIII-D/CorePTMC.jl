## -------- ELEMENT STRUCTURE AND CONSTRUCTOR -------- ##

# Defines structure with flying particle properties (could be an Element, molecule, proton and so on...)

struct Element
    Z::Float64      # charge state
    m::Float64      # mass [kg]
    Zmax::Float64   # atomic number
    Es::Float64     # Sublimation Energy [eV]
    symbol::Symbol
end

function Element(element::Symbol, Z::Float64) # this is a symbol :W, this is a string "W"
    element == :W && return Element(Z, m_amu * 183.84, 74.0, 8.68, element) # if element is W -> create structure Element with properties(Z, ...)
    element == :Kr && return Element(Z, m_amu * 83.79, 36.0, NaN, element)
    element == :Cu && return Element(Z, m_amu * 63.546, 29, 3.49, element) # Cu binding energy from https://doi.org/10.1016/j.vacuum.2011.10.018
    element == :Cr && return Element(Z, m_amu * 52.0, 24.0, 7.41, element)
    element == :Fe && return Element(Z, m_amu * 55.85, 22.0, 4.34, element)
    element == :Ar && return Element(Z, m_amu * 39.95, 18.0, NaN, element)
    element == :Si && return Element(Z, m_amu * 28.09, 14.0, 4.72, element) # Es from RustBCA https://github.com/lcpp-org/RustBCA/blob/main/scripts/materials.py
    element == :Al && return Element(Z, m_amu * 26.98, 13.0, NaN, element) 
    element == :Ne && return Element(Z, m_amu * 20.18, 10.0, NaN, element) 
    element == :C && return Element(Z, m_amu * 12.0, 6.0, 4.12, element)
    element == :B && return Element(Z, m_amu * 10.81, 5.0, NaN, element) 
    element == :Be && return Element(Z, m_amu * 9.01, 4.0, NaN, element)
    element == :He && return Element(Z, m_amu * 4.0, 2.0, NaN, element)
    element == :T && return Element(Z, m_amu * 3.0, 1.0, NaN, element)
    element == :D && return Element(Z, m_amu * 2.0, 1.0, NaN, element)
    element == :H && return Element(Z, m_amu * 1.0, 1.0, NaN, element)

    return error("element $element not available... Need to implement FusionSpecies...")
end

function Element(element::Symbol; Z::Float64 = 0.0) # this is a symbol :W, this is a string "W"
    element == :W && return Element(Z, m_amu * 183.84, 74.0, 8.68, element) # if element is W -> create structure Element with properties(Z, ...)
    element == :Kr && return Element(Z, m_amu * 83.79, 36.0, NaN, element)
    element == :Cu && return Element(Z, m_amu * 63.546, 29, 3.49, element) # Cu binding energy from https://doi.org/10.1016/j.vacuum.2011.10.018
    element == :Cr && return Element(Z, m_amu * 52.0, 24.0, 7.41, element)
    element == :Fe && return Element(Z, m_amu * 55.85, 22.0, 4.34, element)
    element == :Ar && return Element(Z, m_amu * 39.95, 18.0, NaN, element)
    element == :Si && return Element(Z, m_amu * 28.09, 14.0, 4.72, element)
    element == :Al && return Element(Z, m_amu * 26.98, 13.0, NaN, element) 
    element == :Ne && return Element(Z, m_amu * 20.18, 10.0, NaN, element) 
    element == :C && return Element(Z, m_amu * 12.0, 6.0, 4.12, element)
    element == :B && return Element(Z, m_amu * 10.81, 5.0, NaN, element) 
    element == :Be && return Element(Z, m_amu * 9.01, 4.0, NaN, element)
    element == :He && return Element(Z, m_amu * 4.0, 2.0, NaN, element)
    element == :T && return Element(Z, m_amu * 3.0, 1.0, NaN, element)
    element == :D && return Element(Z, m_amu * 2.0, 1.0, NaN, element)
    element == :H && return Element(Z, m_amu * 1.0, 1.0, NaN, element)
    return error("element $element not available... Need to implement FusionSpecies...")
end

get_element(args...; kw...) = Element(args...; kw...) # renaming structure