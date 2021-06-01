#======================================================================
#
# Infix Notation to Reversed Poland Notation Parser
#
#======================================================================

# constants for operator
LEFT_ASSOC = 0
RIGHT_ASSOC = 1

# operators
OPERATORS = {
    '+'    : (0, LEFT_ASSOC),
    '-'    : (0, LEFT_ASSOC),
    '*'    : (5, LEFT_ASSOC),
    '/'    : (5, LEFT_ASSOC),
    '%'    : (5, LEFT_ASSOC),
    '^'    : (10, RIGHT_ASSOC),
    'exp'  : (10, RIGHT_ASSOC),
    'sqrt' : (10, RIGHT_ASSOC)
}

# terms
TERM = {
    '+'    : (1, 1),
    '-'    : (1, 2),
    '*'    : (1, 3),
    '/'    : (1, 4),
    '%'    : (1, 5),
    '^'    : (1, 6),
    'exp'  : (1, 7),
    'sqrt' : (1, 8),
    'Tn'   : (2, 1),
    'Ti'   : (2, 2),
    'Te'   : (2, 3), 
    'M'    : (3, 1)
}

# split characters into tokens
def str_to_infix_with_space(eq):

    # insert ' ' among tokens
    #eq = eq.replace(' ','')
    eq = '('+eq+')'
    eq = eq.replace('T(electron)', 'Te')
    eq = eq.replace('T(ion)', 'Ti')
    eq = eq.replace('T(neutral)', 'Tn')
    eq = eq.replace('e-', 'eminus')
    eq = eq.replace('e+', 'eplus')
    eq = eq.replace('E-', 'eminus')
    eq = eq.replace('E+', 'eplus')
    eq = eq.replace('+', ' + ')
    eq = eq.replace('-', ' - ')
    eq = eq.replace('( + ' , '(+')
    eq = eq.replace('( - ' , '(-')
    eq = eq.replace('**', '^')
    eq = eq.replace('*', ' * ')
    eq = eq.replace('/', ' / ')
    eq = eq.replace('^', ' ^ ')
    eq = eq.replace('exp', ' exp ')
    eq = eq.replace('sqrt', ' sqrt ')
    eq = eq.replace('EXP', ' exp ')
    eq = eq.replace('SQRT', ' sqrt ')
    eq = eq.replace('(', ' ( ')
    eq = eq.replace(')', ' ) ')
    eq = eq.replace('eminus', 'e-')
    eq = eq.replace('eplus', 'e+')

    # split equation by ' '
    eq = eq.split(' ')

    out = []
    for i in range(len(eq)):
        eq[i] = eq[i].lstrip()
        eq[i] = eq[i].rstrip()
        if eq[i] != '':
            out.append(eq[i])
    return out

# operator check
def is_operator(token):
    return token in OPERATORS.keys()

# term check
def is_term(token):
    return token in TERM.keys()

# associativity check
def is_associative(token, assoc):
    if not is_operator(token):
        raise ValueError('Invalid token: %s' % token)
    return OPERATORS[token][1] == assoc

# priority check
def cmp_precedence(token1, token2):
    if not is_operator(token1) or not is_operator(token2):
        raise ValueError('Invalid tokens: %s %s' % (token1, token2))
    return OPERATORS[token1][0] - OPERATORS[token2][0]

# translate infix into RPN
def infix_to_RPN(tokens):
    out = []
    stack = []
    for token in tokens:
        if is_operator(token):
            while len(stack) != 0 and is_operator(stack[-1]):
                if (is_associative(token, LEFT_ASSOC) and
                    cmp_precedence(token, stack[-1]) <= 0) or(
                    is_associative(token, RIGHT_ASSOC) and
                    cmp_precedence(token, stack[-1]) < 0):
                    out.append(stack.pop())
                    continue
                break
            stack.append(token)
        elif token == '(':
            stack.append(token)
        elif token == ')':
            while len(stack) != 0 and stack[-1] != '(':
                out.append(stack.pop())
            stack.pop()
        else:
            out.append(token)
    while len(stack) != 0:
        out.append(stack.pop())
    return out

# for F90 input
def RPN_for_F90(token):
    out = [0]
    label = [0]
    for i in range(len(token)):
        if not is_term(token[i]):
            out.append(token[i])
            label.append(0)
        else:
            out.append(TERM[token[i]][1])
            label.append(TERM[token[i]][0])
    out[0]   = len(out)-1
    label[0] = len(label)-1

    return out,label
