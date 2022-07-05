"""
        Nishank Kankas (EE20B093)
        EE2703 Applied Programming Lab - 2022
        Assignment 2
"""
import sys
import numpy as np

# command to tell starting and ending of circuit
CIRCUIT = '.circuit'
END = '.end'
AC = '.ac'
ground = "GND"

PIE = np.pi     # value of pie

# Is ckt ac or not
ac_circuit = False
circuit_frequency = None


class Comp:
    """
    This class which will be used for components of ckt.

    Note: default frequency is taken as 1e-20 assuming its not ac ckt and hence at this frequency component will behave
    like its dc
    """
    def __init__(self, name=None, n1=None, n2=None, n3=None, n4=None, val=None, is_ac=False, frequency=1e-20,
                 phase=None):
        self.name = name
        self.n1 = n1 if not n1 or n1.isalnum() else sys.exit("Incorrect node name")
        self.n2 = n2 if not n2 or n2.isalnum() else sys.exit("Incorrect node name")
        self.n3 = n3 if not n3 or n3.isalnum() else sys.exit("Incorrect node name")
        self.n4 = n4 if not n4 or n4.isalnum() else sys.exit("Incorrect node name")
        self.value = unit_converter(val)
        self.is_ac = is_ac
        self.frequency = frequency
        self.phase = phase


def unit_converter(value):
    """
    This function just convert the mathematical symbols to the number.
    example:    1k --> 1000     42n --> 0.000000042         2.3G --> 2300000000

    Note: It also convert number in string datatype to float
    """
    units = {
            'y': 1e-24,     'z': 1e-21,     'a': 1e-18,         # yocto       zepto       atto
            'f': 1e-15,     'p': 1e-12,     'n': 1e-9,          # femto       pico        nano
            'u': 1e-6,      'm': 1e-3,      'c': 1e-2,          # micro       milli        centi
            'd': 1e-1,      'k': 1e3,       'M': 1e6,           # deci        kilo        mega
            'G': 1e9,       'T': 1e12,      'P': 1e15,          # giga        tera        peta
            'E': 1e18,      'Z': 1e21,      'Y': 1e24,          # exa         zetta       yotta
           }

    # seeing if there is any unit symbol in end of num
    if value[-1] in units:
        return float(value[:-1]) * units[value-1]
    try:

        # if no, then its just a simple number
        return float(value)
    except ValueError:
        sys.exit("Incorrect unit or value of components")


def open_file():
    """
    This function check if number of argument in command line is correct or not and then it return the file whose name
    is passed in commandline
    """

    # checking if number of argument passes by user is correct or not.
    # There should be two. 1st: Name of python script    2nd: circuit file name
    if len(sys.argv) != 2:
        sys.exit("ERROR: Incorrect numbers of arguments")

    # opening the file and returning it
    try:
        return open(sys.argv[1])
    except FileNotFoundError:
        sys.exit("ERROR!! File not found.")


def get_circuit(content):
    """
    :return: Search for the circuit (if any) in content and return it.
    """
    global ac_circuit

    # initiating parameter containing info about circuit
    circuit_start = circuit_end = False
    start_index = end_index = None
    ac_index = None     # just to remove an warning

    # loop to traverse through content
    for line_index, line in enumerate(content):
        tokens = line.split()

        # Finding start of circuit
        if tokens and tokens[0] == CIRCUIT:
            # Giving error if we had already found circuit start command before
            if circuit_start is True:
                sys.exit(f"ERROR: More then one {CIRCUIT} found")
            else:
                # checking if there is any unwanted stuff after circuit start command
                if len(tokens) > 1 and tokens[1][0] != '#':
                    sys.exit(f"ERROR: Invalid text after {CIRCUIT} in line {line_index + 1}")
                circuit_start = True
                start_index = line_index

        # Finding end of circuit
        if tokens and tokens[0] == END:
            # Giving error if we had already found circuit end command before
            if circuit_end is True:
                sys.exit(f"ERROR: More then one {END} found")
            else:
                # checking if there is any unwanted stuff after circuit end command
                if len(tokens) > 1 and tokens[1][0] != '#':
                    sys.exit(f"ERROR: Invalid text after {END} in line {line_index + 1}")
                circuit_end = True
                end_index = line_index

        # Finding AC command
        if tokens and tokens[0] == AC:
            if ac_circuit is True:
                # same reason what we did before
                sys.exit(f"ERROR: More then one {AC} found")
            elif not circuit_start and not circuit_end:
                sys.exit(f"ERROR: {AC} came before circuit block end")
            else:
                ac_circuit = True
                ac_index = line_index

    # Exiting program if no circuit fond or if circuit block is not valid
    if circuit_start is False:
        sys.exit("Circuit start not Found")

    if circuit_end is False:
        sys.exit("ERROR: Circuit start found but end don't exist")

    if start_index > end_index:
        sys.exit(f"ERROR: Invalid circuit: {CIRCUIT} came after {END}")

    # returning the circuit portion
    return content[start_index+1:end_index] if not ac_circuit else content[start_index+1:end_index]+[content[ac_index]]


def get_components(circuit):
    """
    This function take the circuit as parameter and analysis the branches before returning them
    """
    global circuit_frequency

    if ac_circuit:
        tokens = circuit.pop().split()
        if len(tokens) == 3 or (len(tokens) > 3 and tokens[3][0] == '#'):
            circuit_frequency = 2 * PIE * unit_converter(tokens[2])
        else:
            sys.exit(f"Not valid format of {AC}")

    # initiating dic which will store all components
    components = {
        'R': [],        'I': [],        'L': [],
        'V': [],        'C': [],        'E': [],
        'G': [],        'H': [],        'F': []
    }

    # traversing branches in circuit
    for branch in circuit:
        tokens = branch.split()

        # Will do later
        if tokens and tokens[0][0] in ['E', 'G', 'H', 'F']:
            sys.exit("This version don't support dependent source")

        # if branch contain 4 info
        elif len(tokens) == 4 or (len(tokens) > 4 and tokens[4][0] == '#'):
            if tokens[0][0] == 'R':
                components['R'].append(Comp(name=tokens[0], n1=tokens[1], n2=tokens[2], val=tokens[3]))
            elif tokens[0][0] == 'L':
                components['L'].append(Comp(name=tokens[0], n1=tokens[1], n2=tokens[2], val=tokens[3]))
            elif tokens[0][0] == 'C':
                components['C'].append(Comp(name=tokens[0], n1=tokens[1], n2=tokens[2], val=tokens[3]))
            else:
                sys.exit("incorrect number of token or wrong component name")

        # if branch contain 5 info
        elif len(tokens) == 5 or (len(tokens) > 5 and tokens[5][0] == '#'):
            if tokens[3] == "dc":
                if circuit_frequency:
                    sys.exit("multiple frequency not supported this week")

                if tokens[0][0] == 'V':
                    components['V'].append(Comp(name=tokens[0], n1=tokens[1], n2=tokens[2], val=tokens[4]))
                elif tokens[0][0] == 'I':
                    components['I'].append(Comp(name=tokens[0], n1=tokens[1], n2=tokens[2], val=tokens[4]))
            else:
                sys.exit("Wrong Tokens")

        # if branch contain 6 info
        elif len(tokens) == 6 or (len(tokens) > 6 and tokens[6][0] == '#'):
            if tokens[3] == "ac":
                if not circuit_frequency:
                    sys.exit("ERROR: AC elements given but Frequency not given")

                if tokens[0][0] == 'V':
                    components['V'].append(Comp(name=tokens[0], n1=tokens[1], n2=tokens[2], val=tokens[4], is_ac=True,
                                                frequency=circuit_frequency, phase=tokens[5]))
                elif tokens[0][0] == 'I':
                    components['I'].append(Comp(name=tokens[0], n1=tokens[1], n2=tokens[2], val=tokens[4], is_ac=True,
                                                frequency=circuit_frequency, phase=tokens[5]))
            else:
                sys.exit("Wrong tokens")

        elif not tokens:
            sys.exit("ERROR: An empty line found in circuit")

        # if its neither of above, them the line must be comment or it is an invalid line
        elif tokens[0][0] != '#':
            sys.exit("ERROR: Invalid line found in circuit")

    return components


def correct_val(components):
    if not ac_circuit:
        for elem in components["L"]:
            elem.value = 0
        for elem in components["C"]:
            elem.value = np.inf

    elif ac_circuit:
        for elem in components["V"]:
            elem.value = elem.value/2*np.exp(1j*float(elem.phase)*np.pi/180)

        for elem in components["I"]:
            elem.value = elem.value/2*np.exp(1j*elem.phase*PIE/180)

        for elem in components["C"]:
            elem.value = 1/(1j * elem.value * circuit_frequency)

        for elem in components["L"]:
            elem.value = 1j * elem.value * circuit_frequency

    return components


def get_list(components):
    n = []
    v = []

    for comp in components['V']:
        v.append(comp.name)

    for comp_type in components:
        for comp in components[comp_type]:
            if comp.n1 and comp.n1 not in n:
                n.append(comp.n1)
            if comp.n2 and comp.n2 not in n:
                n.append(comp.n2)
            if comp.n3 and comp.n3 not in n:
                n.append(comp.n3)
            if comp.n4 and comp.n4 not in n:
                n.append(comp.n4)
    return n, v


def get_matrix(node_list, volt_list, components):
    m_list = node_list + volt_list
    m_size = len(m_list)
    m = np.zeros((m_size, m_size), dtype='complex')
    b = np.zeros((m_size, 1), dtype='complex')

    for comp_type in ['R', 'C', 'L']:
        for comp in components[comp_type]:
            gm = 1/comp.value
            m_i = m_list.index(comp.n1)
            m_j = m_list.index(comp.n2)
            m[m_i, m_i] = m[m_i, m_i] + gm
            m[m_j, m_j] = m[m_j, m_j] + gm
            m[m_i, m_j] = m[m_i, m_j] - gm
            m[m_j, m_i] = m[m_j, m_i] - gm

    for comp in components['I']:
        m_i = m_list.index(comp.n1)
        m_j = m_list.index(comp.n2)
        b[m_i] = b[m_i] + comp.value
        b[m_j] = b[m_j] - comp.value

    for comp in components['V']:
        m_i = m_list.index(comp.n1)
        m_j = m_list.index(comp.n2)
        b_i = m_list.index(comp.name)
        m[m_i, b_i] = m[m_i, b_i] + 1
        m[b_i, m_i] = m[b_i, m_i] + 1
        m[m_j, b_i] = m[m_j, b_i] - 1
        m[b_i, m_j] = m[b_i, m_j] - 1
        b[b_i] = b[b_i] + comp.value

    global ground
    if "GND" not in node_list:
        ground = node_list[0]
        print(node_list[0], "is chosen as ground")

    round_index = m_list.index(ground)
    m = np.delete(m, round_index, 0)
    m = np.delete(m, round_index, 1)
    b = np.delete(b, round_index, 0)
    m_list.remove(ground)
    node_list.remove(ground)
    return m, b


def print_result(b, x, node_list, components):
    print("Voltage at node ", ground, " is ", 0)
    for i in range(len(b)):
        if i < len(node_list):
            print(f"Voltage at node {node_list[i]} is {x[i]}")
        elif i < len(node_list) + len(components['I']):
            print(f"Current through {components['I'][i - len(node_list)].name} is {x[i]}")
        else:
            print(f"Current through {components['V'][i - len(node_list) - len(components['I'])].name} is {x[i]}")


def main():
    file = open_file()                      # opening file
    content = file.readlines()              # reading content in file
    circuit = get_circuit(content)          # Extracting the circuit (note: get_circuit() will return in reverse order)
    file.close()                            # Closing the file
    components = get_components(circuit)    # Extracting components in circuit
    components = correct_val(components)    # Printing branches in reverse order to print whole circuit in req format
    node_list, volt_list = get_list(components)              # getting list of nodes and current
    m, b = get_matrix(node_list, volt_list, components)      # finding m and b matrix
    x = np.linalg.solve(m, b)                                   # solving it using numpy linear algebra solver
    print_result(b, x, node_list, components)                   # Finally print the result


if __name__ == '__main__':
    main()
