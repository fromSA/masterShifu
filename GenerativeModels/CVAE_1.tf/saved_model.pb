��
��
:
Add
x"T
y"T
z"T"
Ttype:
2	
B
AddV2
x"T
y"T
z"T"
Ttype:
2	��
B
AssignVariableOp
resource
value"dtype"
dtypetype�
~
BiasAdd

value"T	
bias"T
output"T" 
Ttype:
2	"-
data_formatstringNHWC:
NHWCNCHW
h
ConcatV2
values"T*N
axis"Tidx
output"T"
Nint(0"	
Ttype"
Tidxtype0:
2	
8
Const
output"dtype"
valuetensor"
dtypetype
,
Exp
x"T
y"T"
Ttype:

2
.
Identity

input"T
output"T"	
Ttype
q
MatMul
a"T
b"T
product"T"
transpose_abool( "
transpose_bbool( "
Ttype:

2	
e
MergeV2Checkpoints
checkpoint_prefixes
destination_prefix"
delete_old_dirsbool(�
=
Mul
x"T
y"T
z"T"
Ttype:
2	�

NoOp
M
Pack
values"T*N
output"T"
Nint(0"	
Ttype"
axisint 
C
Placeholder
output"dtype"
dtypetype"
shapeshape:
�
RandomStandardNormal

shape"T
output"dtype"
seedint "
seed2int "
dtypetype:
2"
Ttype:
2	�
@
ReadVariableOp
resource
value"dtype"
dtypetype�
E
Relu
features"T
activations"T"
Ttype:
2	
o
	RestoreV2

prefix
tensor_names
shape_and_slices
tensors2dtypes"
dtypes
list(type)(0�
l
SaveV2

prefix
tensor_names
shape_and_slices
tensors2dtypes"
dtypes
list(type)(0�
?
Select
	condition

t"T
e"T
output"T"	
Ttype
P
Shape

input"T
output"out_type"	
Ttype"
out_typetype0:
2	
H
ShardedFilename
basename	
shard

num_shards
filename
�
StatefulPartitionedCall
args2Tin
output2Tout"
Tin
list(type)("
Tout
list(type)("	
ffunc"
configstring "
config_protostring "
executor_typestring �
@
StaticRegexFullMatch	
input

output
"
patternstring
�
StridedSlice

input"T
begin"Index
end"Index
strides"Index
output"T"	
Ttype"
Indextype:
2	"

begin_maskint "
end_maskint "
ellipsis_maskint "
new_axis_maskint "
shrink_axis_maskint 
N

StringJoin
inputs*N

output"
Nint(0"
	separatorstring 
�
VarHandleOp
resource"
	containerstring "
shared_namestring "
dtypetype"
shapeshape"#
allowed_deviceslist(string)
 �"serve*2.4.02v2.4.0-rc4-71-g582c8d236cb8��
f
	Adam/iterVarHandleOp*
_output_shapes
: *
dtype0	*
shape: *
shared_name	Adam/iter
_
Adam/iter/Read/ReadVariableOpReadVariableOp	Adam/iter*
_output_shapes
: *
dtype0	
j
Adam/beta_1VarHandleOp*
_output_shapes
: *
dtype0*
shape: *
shared_nameAdam/beta_1
c
Adam/beta_1/Read/ReadVariableOpReadVariableOpAdam/beta_1*
_output_shapes
: *
dtype0
j
Adam/beta_2VarHandleOp*
_output_shapes
: *
dtype0*
shape: *
shared_nameAdam/beta_2
c
Adam/beta_2/Read/ReadVariableOpReadVariableOpAdam/beta_2*
_output_shapes
: *
dtype0
h

Adam/decayVarHandleOp*
_output_shapes
: *
dtype0*
shape: *
shared_name
Adam/decay
a
Adam/decay/Read/ReadVariableOpReadVariableOp
Adam/decay*
_output_shapes
: *
dtype0
x
Adam/learning_rateVarHandleOp*
_output_shapes
: *
dtype0*
shape: *#
shared_nameAdam/learning_rate
q
&Adam/learning_rate/Read/ReadVariableOpReadVariableOpAdam/learning_rate*
_output_shapes
: *
dtype0
~
encoder_l2/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:
*"
shared_nameencoder_l2/kernel
w
%encoder_l2/kernel/Read/ReadVariableOpReadVariableOpencoder_l2/kernel*
_output_shapes

:
*
dtype0
v
encoder_l2/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:
* 
shared_nameencoder_l2/bias
o
#encoder_l2/bias/Read/ReadVariableOpReadVariableOpencoder_l2/bias*
_output_shapes
:
*
dtype0
~
encoder_l3/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:
*"
shared_nameencoder_l3/kernel
w
%encoder_l3/kernel/Read/ReadVariableOpReadVariableOpencoder_l3/kernel*
_output_shapes

:
*
dtype0
v
encoder_l3/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:* 
shared_nameencoder_l3/bias
o
#encoder_l3/bias/Read/ReadVariableOpReadVariableOpencoder_l3/bias*
_output_shapes
:*
dtype0
~
encoder_l4/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:*"
shared_nameencoder_l4/kernel
w
%encoder_l4/kernel/Read/ReadVariableOpReadVariableOpencoder_l4/kernel*
_output_shapes

:*
dtype0
v
encoder_l4/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:* 
shared_nameencoder_l4/bias
o
#encoder_l4/bias/Read/ReadVariableOpReadVariableOpencoder_l4/bias*
_output_shapes
:*
dtype0
v
z_mean/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:*
shared_namez_mean/kernel
o
!z_mean/kernel/Read/ReadVariableOpReadVariableOpz_mean/kernel*
_output_shapes

:*
dtype0
n
z_mean/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:*
shared_namez_mean/bias
g
z_mean/bias/Read/ReadVariableOpReadVariableOpz_mean/bias*
_output_shapes
:*
dtype0
|
z_log_var/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:*!
shared_namez_log_var/kernel
u
$z_log_var/kernel/Read/ReadVariableOpReadVariableOpz_log_var/kernel*
_output_shapes

:*
dtype0
t
z_log_var/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:*
shared_namez_log_var/bias
m
"z_log_var/bias/Read/ReadVariableOpReadVariableOpz_log_var/bias*
_output_shapes
:*
dtype0
~
decoder_l2/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:*"
shared_namedecoder_l2/kernel
w
%decoder_l2/kernel/Read/ReadVariableOpReadVariableOpdecoder_l2/kernel*
_output_shapes

:*
dtype0
v
decoder_l2/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:* 
shared_namedecoder_l2/bias
o
#decoder_l2/bias/Read/ReadVariableOpReadVariableOpdecoder_l2/bias*
_output_shapes
:*
dtype0
~
decoder_l3/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:*"
shared_namedecoder_l3/kernel
w
%decoder_l3/kernel/Read/ReadVariableOpReadVariableOpdecoder_l3/kernel*
_output_shapes

:*
dtype0
v
decoder_l3/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:* 
shared_namedecoder_l3/bias
o
#decoder_l3/bias/Read/ReadVariableOpReadVariableOpdecoder_l3/bias*
_output_shapes
:*
dtype0
~
decoder_l4/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:
*"
shared_namedecoder_l4/kernel
w
%decoder_l4/kernel/Read/ReadVariableOpReadVariableOpdecoder_l4/kernel*
_output_shapes

:
*
dtype0
v
decoder_l4/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:
* 
shared_namedecoder_l4/bias
o
#decoder_l4/bias/Read/ReadVariableOpReadVariableOpdecoder_l4/bias*
_output_shapes
:
*
dtype0
�
output_layer/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:
*$
shared_nameoutput_layer/kernel
{
'output_layer/kernel/Read/ReadVariableOpReadVariableOpoutput_layer/kernel*
_output_shapes

:
*
dtype0
z
output_layer/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:*"
shared_nameoutput_layer/bias
s
%output_layer/bias/Read/ReadVariableOpReadVariableOpoutput_layer/bias*
_output_shapes
:*
dtype0
^
totalVarHandleOp*
_output_shapes
: *
dtype0*
shape: *
shared_nametotal
W
total/Read/ReadVariableOpReadVariableOptotal*
_output_shapes
: *
dtype0
^
countVarHandleOp*
_output_shapes
: *
dtype0*
shape: *
shared_namecount
W
count/Read/ReadVariableOpReadVariableOpcount*
_output_shapes
: *
dtype0
�
Adam/encoder_l2/kernel/mVarHandleOp*
_output_shapes
: *
dtype0*
shape
:
*)
shared_nameAdam/encoder_l2/kernel/m
�
,Adam/encoder_l2/kernel/m/Read/ReadVariableOpReadVariableOpAdam/encoder_l2/kernel/m*
_output_shapes

:
*
dtype0
�
Adam/encoder_l2/bias/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:
*'
shared_nameAdam/encoder_l2/bias/m
}
*Adam/encoder_l2/bias/m/Read/ReadVariableOpReadVariableOpAdam/encoder_l2/bias/m*
_output_shapes
:
*
dtype0
�
Adam/encoder_l3/kernel/mVarHandleOp*
_output_shapes
: *
dtype0*
shape
:
*)
shared_nameAdam/encoder_l3/kernel/m
�
,Adam/encoder_l3/kernel/m/Read/ReadVariableOpReadVariableOpAdam/encoder_l3/kernel/m*
_output_shapes

:
*
dtype0
�
Adam/encoder_l3/bias/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:*'
shared_nameAdam/encoder_l3/bias/m
}
*Adam/encoder_l3/bias/m/Read/ReadVariableOpReadVariableOpAdam/encoder_l3/bias/m*
_output_shapes
:*
dtype0
�
Adam/encoder_l4/kernel/mVarHandleOp*
_output_shapes
: *
dtype0*
shape
:*)
shared_nameAdam/encoder_l4/kernel/m
�
,Adam/encoder_l4/kernel/m/Read/ReadVariableOpReadVariableOpAdam/encoder_l4/kernel/m*
_output_shapes

:*
dtype0
�
Adam/encoder_l4/bias/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:*'
shared_nameAdam/encoder_l4/bias/m
}
*Adam/encoder_l4/bias/m/Read/ReadVariableOpReadVariableOpAdam/encoder_l4/bias/m*
_output_shapes
:*
dtype0
�
Adam/z_mean/kernel/mVarHandleOp*
_output_shapes
: *
dtype0*
shape
:*%
shared_nameAdam/z_mean/kernel/m
}
(Adam/z_mean/kernel/m/Read/ReadVariableOpReadVariableOpAdam/z_mean/kernel/m*
_output_shapes

:*
dtype0
|
Adam/z_mean/bias/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:*#
shared_nameAdam/z_mean/bias/m
u
&Adam/z_mean/bias/m/Read/ReadVariableOpReadVariableOpAdam/z_mean/bias/m*
_output_shapes
:*
dtype0
�
Adam/z_log_var/kernel/mVarHandleOp*
_output_shapes
: *
dtype0*
shape
:*(
shared_nameAdam/z_log_var/kernel/m
�
+Adam/z_log_var/kernel/m/Read/ReadVariableOpReadVariableOpAdam/z_log_var/kernel/m*
_output_shapes

:*
dtype0
�
Adam/z_log_var/bias/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:*&
shared_nameAdam/z_log_var/bias/m
{
)Adam/z_log_var/bias/m/Read/ReadVariableOpReadVariableOpAdam/z_log_var/bias/m*
_output_shapes
:*
dtype0
�
Adam/decoder_l2/kernel/mVarHandleOp*
_output_shapes
: *
dtype0*
shape
:*)
shared_nameAdam/decoder_l2/kernel/m
�
,Adam/decoder_l2/kernel/m/Read/ReadVariableOpReadVariableOpAdam/decoder_l2/kernel/m*
_output_shapes

:*
dtype0
�
Adam/decoder_l2/bias/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:*'
shared_nameAdam/decoder_l2/bias/m
}
*Adam/decoder_l2/bias/m/Read/ReadVariableOpReadVariableOpAdam/decoder_l2/bias/m*
_output_shapes
:*
dtype0
�
Adam/decoder_l3/kernel/mVarHandleOp*
_output_shapes
: *
dtype0*
shape
:*)
shared_nameAdam/decoder_l3/kernel/m
�
,Adam/decoder_l3/kernel/m/Read/ReadVariableOpReadVariableOpAdam/decoder_l3/kernel/m*
_output_shapes

:*
dtype0
�
Adam/decoder_l3/bias/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:*'
shared_nameAdam/decoder_l3/bias/m
}
*Adam/decoder_l3/bias/m/Read/ReadVariableOpReadVariableOpAdam/decoder_l3/bias/m*
_output_shapes
:*
dtype0
�
Adam/decoder_l4/kernel/mVarHandleOp*
_output_shapes
: *
dtype0*
shape
:
*)
shared_nameAdam/decoder_l4/kernel/m
�
,Adam/decoder_l4/kernel/m/Read/ReadVariableOpReadVariableOpAdam/decoder_l4/kernel/m*
_output_shapes

:
*
dtype0
�
Adam/decoder_l4/bias/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:
*'
shared_nameAdam/decoder_l4/bias/m
}
*Adam/decoder_l4/bias/m/Read/ReadVariableOpReadVariableOpAdam/decoder_l4/bias/m*
_output_shapes
:
*
dtype0
�
Adam/output_layer/kernel/mVarHandleOp*
_output_shapes
: *
dtype0*
shape
:
*+
shared_nameAdam/output_layer/kernel/m
�
.Adam/output_layer/kernel/m/Read/ReadVariableOpReadVariableOpAdam/output_layer/kernel/m*
_output_shapes

:
*
dtype0
�
Adam/output_layer/bias/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:*)
shared_nameAdam/output_layer/bias/m
�
,Adam/output_layer/bias/m/Read/ReadVariableOpReadVariableOpAdam/output_layer/bias/m*
_output_shapes
:*
dtype0
�
Adam/encoder_l2/kernel/vVarHandleOp*
_output_shapes
: *
dtype0*
shape
:
*)
shared_nameAdam/encoder_l2/kernel/v
�
,Adam/encoder_l2/kernel/v/Read/ReadVariableOpReadVariableOpAdam/encoder_l2/kernel/v*
_output_shapes

:
*
dtype0
�
Adam/encoder_l2/bias/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:
*'
shared_nameAdam/encoder_l2/bias/v
}
*Adam/encoder_l2/bias/v/Read/ReadVariableOpReadVariableOpAdam/encoder_l2/bias/v*
_output_shapes
:
*
dtype0
�
Adam/encoder_l3/kernel/vVarHandleOp*
_output_shapes
: *
dtype0*
shape
:
*)
shared_nameAdam/encoder_l3/kernel/v
�
,Adam/encoder_l3/kernel/v/Read/ReadVariableOpReadVariableOpAdam/encoder_l3/kernel/v*
_output_shapes

:
*
dtype0
�
Adam/encoder_l3/bias/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:*'
shared_nameAdam/encoder_l3/bias/v
}
*Adam/encoder_l3/bias/v/Read/ReadVariableOpReadVariableOpAdam/encoder_l3/bias/v*
_output_shapes
:*
dtype0
�
Adam/encoder_l4/kernel/vVarHandleOp*
_output_shapes
: *
dtype0*
shape
:*)
shared_nameAdam/encoder_l4/kernel/v
�
,Adam/encoder_l4/kernel/v/Read/ReadVariableOpReadVariableOpAdam/encoder_l4/kernel/v*
_output_shapes

:*
dtype0
�
Adam/encoder_l4/bias/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:*'
shared_nameAdam/encoder_l4/bias/v
}
*Adam/encoder_l4/bias/v/Read/ReadVariableOpReadVariableOpAdam/encoder_l4/bias/v*
_output_shapes
:*
dtype0
�
Adam/z_mean/kernel/vVarHandleOp*
_output_shapes
: *
dtype0*
shape
:*%
shared_nameAdam/z_mean/kernel/v
}
(Adam/z_mean/kernel/v/Read/ReadVariableOpReadVariableOpAdam/z_mean/kernel/v*
_output_shapes

:*
dtype0
|
Adam/z_mean/bias/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:*#
shared_nameAdam/z_mean/bias/v
u
&Adam/z_mean/bias/v/Read/ReadVariableOpReadVariableOpAdam/z_mean/bias/v*
_output_shapes
:*
dtype0
�
Adam/z_log_var/kernel/vVarHandleOp*
_output_shapes
: *
dtype0*
shape
:*(
shared_nameAdam/z_log_var/kernel/v
�
+Adam/z_log_var/kernel/v/Read/ReadVariableOpReadVariableOpAdam/z_log_var/kernel/v*
_output_shapes

:*
dtype0
�
Adam/z_log_var/bias/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:*&
shared_nameAdam/z_log_var/bias/v
{
)Adam/z_log_var/bias/v/Read/ReadVariableOpReadVariableOpAdam/z_log_var/bias/v*
_output_shapes
:*
dtype0
�
Adam/decoder_l2/kernel/vVarHandleOp*
_output_shapes
: *
dtype0*
shape
:*)
shared_nameAdam/decoder_l2/kernel/v
�
,Adam/decoder_l2/kernel/v/Read/ReadVariableOpReadVariableOpAdam/decoder_l2/kernel/v*
_output_shapes

:*
dtype0
�
Adam/decoder_l2/bias/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:*'
shared_nameAdam/decoder_l2/bias/v
}
*Adam/decoder_l2/bias/v/Read/ReadVariableOpReadVariableOpAdam/decoder_l2/bias/v*
_output_shapes
:*
dtype0
�
Adam/decoder_l3/kernel/vVarHandleOp*
_output_shapes
: *
dtype0*
shape
:*)
shared_nameAdam/decoder_l3/kernel/v
�
,Adam/decoder_l3/kernel/v/Read/ReadVariableOpReadVariableOpAdam/decoder_l3/kernel/v*
_output_shapes

:*
dtype0
�
Adam/decoder_l3/bias/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:*'
shared_nameAdam/decoder_l3/bias/v
}
*Adam/decoder_l3/bias/v/Read/ReadVariableOpReadVariableOpAdam/decoder_l3/bias/v*
_output_shapes
:*
dtype0
�
Adam/decoder_l4/kernel/vVarHandleOp*
_output_shapes
: *
dtype0*
shape
:
*)
shared_nameAdam/decoder_l4/kernel/v
�
,Adam/decoder_l4/kernel/v/Read/ReadVariableOpReadVariableOpAdam/decoder_l4/kernel/v*
_output_shapes

:
*
dtype0
�
Adam/decoder_l4/bias/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:
*'
shared_nameAdam/decoder_l4/bias/v
}
*Adam/decoder_l4/bias/v/Read/ReadVariableOpReadVariableOpAdam/decoder_l4/bias/v*
_output_shapes
:
*
dtype0
�
Adam/output_layer/kernel/vVarHandleOp*
_output_shapes
: *
dtype0*
shape
:
*+
shared_nameAdam/output_layer/kernel/v
�
.Adam/output_layer/kernel/v/Read/ReadVariableOpReadVariableOpAdam/output_layer/kernel/v*
_output_shapes

:
*
dtype0
�
Adam/output_layer/bias/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:*)
shared_nameAdam/output_layer/bias/v
�
,Adam/output_layer/bias/v/Read/ReadVariableOpReadVariableOpAdam/output_layer/bias/v*
_output_shapes
:*
dtype0

NoOpNoOp
�\
ConstConst"/device:CPU:0*
_output_shapes
: *
dtype0*�[
value�[B�[ B�[
�
encoder
decoder
	optimizer
loss
	variables
trainable_variables
regularization_losses
	keras_api
	
signatures
�

layer-0
layer_with_weights-0
layer-1
layer_with_weights-1
layer-2
layer_with_weights-2
layer-3
layer_with_weights-3
layer-4
layer_with_weights-4
layer-5
layer-6
	variables
trainable_variables
regularization_losses
	keras_api
�
layer-0
layer_with_weights-0
layer-1
layer_with_weights-1
layer-2
layer_with_weights-2
layer-3
layer_with_weights-3
layer-4
	variables
trainable_variables
regularization_losses
	keras_api
�
iter

beta_1

 beta_2
	!decay
"learning_rate#m�$m�%m�&m�'m�(m�)m�*m�+m�,m�-m�.m�/m�0m�1m�2m�3m�4m�#v�$v�%v�&v�'v�(v�)v�*v�+v�,v�-v�.v�/v�0v�1v�2v�3v�4v�
 
�
#0
$1
%2
&3
'4
(5
)6
*7
+8
,9
-10
.11
/12
013
114
215
316
417
�
#0
$1
%2
&3
'4
(5
)6
*7
+8
,9
-10
.11
/12
013
114
215
316
417
 
�
	variables
5metrics
6layer_regularization_losses
trainable_variables
7non_trainable_variables

8layers
9layer_metrics
regularization_losses
 
 
h

#kernel
$bias
:	variables
;trainable_variables
<regularization_losses
=	keras_api
h

%kernel
&bias
>	variables
?trainable_variables
@regularization_losses
A	keras_api
h

'kernel
(bias
B	variables
Ctrainable_variables
Dregularization_losses
E	keras_api
h

)kernel
*bias
F	variables
Gtrainable_variables
Hregularization_losses
I	keras_api
h

+kernel
,bias
J	variables
Ktrainable_variables
Lregularization_losses
M	keras_api
R
N	variables
Otrainable_variables
Pregularization_losses
Q	keras_api
F
#0
$1
%2
&3
'4
(5
)6
*7
+8
,9
F
#0
$1
%2
&3
'4
(5
)6
*7
+8
,9
 
�
	variables
Rmetrics
Slayer_regularization_losses
trainable_variables
Tnon_trainable_variables

Ulayers
Vlayer_metrics
regularization_losses
 
h

-kernel
.bias
W	variables
Xtrainable_variables
Yregularization_losses
Z	keras_api
h

/kernel
0bias
[	variables
\trainable_variables
]regularization_losses
^	keras_api
h

1kernel
2bias
_	variables
`trainable_variables
aregularization_losses
b	keras_api
h

3kernel
4bias
c	variables
dtrainable_variables
eregularization_losses
f	keras_api
8
-0
.1
/2
03
14
25
36
47
8
-0
.1
/2
03
14
25
36
47
 
�
	variables
gmetrics
hlayer_regularization_losses
trainable_variables
inon_trainable_variables

jlayers
klayer_metrics
regularization_losses
HF
VARIABLE_VALUE	Adam/iter)optimizer/iter/.ATTRIBUTES/VARIABLE_VALUE
LJ
VARIABLE_VALUEAdam/beta_1+optimizer/beta_1/.ATTRIBUTES/VARIABLE_VALUE
LJ
VARIABLE_VALUEAdam/beta_2+optimizer/beta_2/.ATTRIBUTES/VARIABLE_VALUE
JH
VARIABLE_VALUE
Adam/decay*optimizer/decay/.ATTRIBUTES/VARIABLE_VALUE
ZX
VARIABLE_VALUEAdam/learning_rate2optimizer/learning_rate/.ATTRIBUTES/VARIABLE_VALUE
MK
VARIABLE_VALUEencoder_l2/kernel&variables/0/.ATTRIBUTES/VARIABLE_VALUE
KI
VARIABLE_VALUEencoder_l2/bias&variables/1/.ATTRIBUTES/VARIABLE_VALUE
MK
VARIABLE_VALUEencoder_l3/kernel&variables/2/.ATTRIBUTES/VARIABLE_VALUE
KI
VARIABLE_VALUEencoder_l3/bias&variables/3/.ATTRIBUTES/VARIABLE_VALUE
MK
VARIABLE_VALUEencoder_l4/kernel&variables/4/.ATTRIBUTES/VARIABLE_VALUE
KI
VARIABLE_VALUEencoder_l4/bias&variables/5/.ATTRIBUTES/VARIABLE_VALUE
IG
VARIABLE_VALUEz_mean/kernel&variables/6/.ATTRIBUTES/VARIABLE_VALUE
GE
VARIABLE_VALUEz_mean/bias&variables/7/.ATTRIBUTES/VARIABLE_VALUE
LJ
VARIABLE_VALUEz_log_var/kernel&variables/8/.ATTRIBUTES/VARIABLE_VALUE
JH
VARIABLE_VALUEz_log_var/bias&variables/9/.ATTRIBUTES/VARIABLE_VALUE
NL
VARIABLE_VALUEdecoder_l2/kernel'variables/10/.ATTRIBUTES/VARIABLE_VALUE
LJ
VARIABLE_VALUEdecoder_l2/bias'variables/11/.ATTRIBUTES/VARIABLE_VALUE
NL
VARIABLE_VALUEdecoder_l3/kernel'variables/12/.ATTRIBUTES/VARIABLE_VALUE
LJ
VARIABLE_VALUEdecoder_l3/bias'variables/13/.ATTRIBUTES/VARIABLE_VALUE
NL
VARIABLE_VALUEdecoder_l4/kernel'variables/14/.ATTRIBUTES/VARIABLE_VALUE
LJ
VARIABLE_VALUEdecoder_l4/bias'variables/15/.ATTRIBUTES/VARIABLE_VALUE
PN
VARIABLE_VALUEoutput_layer/kernel'variables/16/.ATTRIBUTES/VARIABLE_VALUE
NL
VARIABLE_VALUEoutput_layer/bias'variables/17/.ATTRIBUTES/VARIABLE_VALUE

l0
 
 

0
1
 

#0
$1

#0
$1
 
�
:	variables
mmetrics
nlayer_regularization_losses
;trainable_variables
onon_trainable_variables

players
qlayer_metrics
<regularization_losses

%0
&1

%0
&1
 
�
>	variables
rmetrics
slayer_regularization_losses
?trainable_variables
tnon_trainable_variables

ulayers
vlayer_metrics
@regularization_losses

'0
(1

'0
(1
 
�
B	variables
wmetrics
xlayer_regularization_losses
Ctrainable_variables
ynon_trainable_variables

zlayers
{layer_metrics
Dregularization_losses

)0
*1

)0
*1
 
�
F	variables
|metrics
}layer_regularization_losses
Gtrainable_variables
~non_trainable_variables

layers
�layer_metrics
Hregularization_losses

+0
,1

+0
,1
 
�
J	variables
�metrics
 �layer_regularization_losses
Ktrainable_variables
�non_trainable_variables
�layers
�layer_metrics
Lregularization_losses
 
 
 
�
N	variables
�metrics
 �layer_regularization_losses
Otrainable_variables
�non_trainable_variables
�layers
�layer_metrics
Pregularization_losses
 
 
 
1

0
1
2
3
4
5
6
 

-0
.1

-0
.1
 
�
W	variables
�metrics
 �layer_regularization_losses
Xtrainable_variables
�non_trainable_variables
�layers
�layer_metrics
Yregularization_losses

/0
01

/0
01
 
�
[	variables
�metrics
 �layer_regularization_losses
\trainable_variables
�non_trainable_variables
�layers
�layer_metrics
]regularization_losses

10
21

10
21
 
�
_	variables
�metrics
 �layer_regularization_losses
`trainable_variables
�non_trainable_variables
�layers
�layer_metrics
aregularization_losses

30
41

30
41
 
�
c	variables
�metrics
 �layer_regularization_losses
dtrainable_variables
�non_trainable_variables
�layers
�layer_metrics
eregularization_losses
 
 
 
#
0
1
2
3
4
 
8

�total

�count
�	variables
�	keras_api
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
OM
VARIABLE_VALUEtotal4keras_api/metrics/0/total/.ATTRIBUTES/VARIABLE_VALUE
OM
VARIABLE_VALUEcount4keras_api/metrics/0/count/.ATTRIBUTES/VARIABLE_VALUE

�0
�1

�	variables
pn
VARIABLE_VALUEAdam/encoder_l2/kernel/mBvariables/0/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
nl
VARIABLE_VALUEAdam/encoder_l2/bias/mBvariables/1/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
pn
VARIABLE_VALUEAdam/encoder_l3/kernel/mBvariables/2/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
nl
VARIABLE_VALUEAdam/encoder_l3/bias/mBvariables/3/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
pn
VARIABLE_VALUEAdam/encoder_l4/kernel/mBvariables/4/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
nl
VARIABLE_VALUEAdam/encoder_l4/bias/mBvariables/5/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
lj
VARIABLE_VALUEAdam/z_mean/kernel/mBvariables/6/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
jh
VARIABLE_VALUEAdam/z_mean/bias/mBvariables/7/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
om
VARIABLE_VALUEAdam/z_log_var/kernel/mBvariables/8/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
mk
VARIABLE_VALUEAdam/z_log_var/bias/mBvariables/9/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
qo
VARIABLE_VALUEAdam/decoder_l2/kernel/mCvariables/10/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
om
VARIABLE_VALUEAdam/decoder_l2/bias/mCvariables/11/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
qo
VARIABLE_VALUEAdam/decoder_l3/kernel/mCvariables/12/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
om
VARIABLE_VALUEAdam/decoder_l3/bias/mCvariables/13/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
qo
VARIABLE_VALUEAdam/decoder_l4/kernel/mCvariables/14/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
om
VARIABLE_VALUEAdam/decoder_l4/bias/mCvariables/15/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
sq
VARIABLE_VALUEAdam/output_layer/kernel/mCvariables/16/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
qo
VARIABLE_VALUEAdam/output_layer/bias/mCvariables/17/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
pn
VARIABLE_VALUEAdam/encoder_l2/kernel/vBvariables/0/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
nl
VARIABLE_VALUEAdam/encoder_l2/bias/vBvariables/1/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
pn
VARIABLE_VALUEAdam/encoder_l3/kernel/vBvariables/2/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
nl
VARIABLE_VALUEAdam/encoder_l3/bias/vBvariables/3/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
pn
VARIABLE_VALUEAdam/encoder_l4/kernel/vBvariables/4/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
nl
VARIABLE_VALUEAdam/encoder_l4/bias/vBvariables/5/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
lj
VARIABLE_VALUEAdam/z_mean/kernel/vBvariables/6/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
jh
VARIABLE_VALUEAdam/z_mean/bias/vBvariables/7/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
om
VARIABLE_VALUEAdam/z_log_var/kernel/vBvariables/8/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
mk
VARIABLE_VALUEAdam/z_log_var/bias/vBvariables/9/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
qo
VARIABLE_VALUEAdam/decoder_l2/kernel/vCvariables/10/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
om
VARIABLE_VALUEAdam/decoder_l2/bias/vCvariables/11/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
qo
VARIABLE_VALUEAdam/decoder_l3/kernel/vCvariables/12/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
om
VARIABLE_VALUEAdam/decoder_l3/bias/vCvariables/13/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
qo
VARIABLE_VALUEAdam/decoder_l4/kernel/vCvariables/14/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
om
VARIABLE_VALUEAdam/decoder_l4/bias/vCvariables/15/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
sq
VARIABLE_VALUEAdam/output_layer/kernel/vCvariables/16/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
qo
VARIABLE_VALUEAdam/output_layer/bias/vCvariables/17/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
z
serving_default_input_1Placeholder*'
_output_shapes
:���������*
dtype0*
shape:���������
�
StatefulPartitionedCallStatefulPartitionedCallserving_default_input_1encoder_l2/kernelencoder_l2/biasencoder_l3/kernelencoder_l3/biasencoder_l4/kernelencoder_l4/biasz_mean/kernelz_mean/biasz_log_var/kernelz_log_var/biasdecoder_l2/kerneldecoder_l2/biasdecoder_l3/kerneldecoder_l3/biasdecoder_l4/kerneldecoder_l4/biasoutput_layer/kerneloutput_layer/bias*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*4
_read_only_resource_inputs
	
*-
config_proto

CPU

GPU 2J 8� *-
f(R&
$__inference_signature_wrapper_142095
O
saver_filenamePlaceholder*
_output_shapes
: *
dtype0*
shape: 
�
StatefulPartitionedCall_1StatefulPartitionedCallsaver_filenameAdam/iter/Read/ReadVariableOpAdam/beta_1/Read/ReadVariableOpAdam/beta_2/Read/ReadVariableOpAdam/decay/Read/ReadVariableOp&Adam/learning_rate/Read/ReadVariableOp%encoder_l2/kernel/Read/ReadVariableOp#encoder_l2/bias/Read/ReadVariableOp%encoder_l3/kernel/Read/ReadVariableOp#encoder_l3/bias/Read/ReadVariableOp%encoder_l4/kernel/Read/ReadVariableOp#encoder_l4/bias/Read/ReadVariableOp!z_mean/kernel/Read/ReadVariableOpz_mean/bias/Read/ReadVariableOp$z_log_var/kernel/Read/ReadVariableOp"z_log_var/bias/Read/ReadVariableOp%decoder_l2/kernel/Read/ReadVariableOp#decoder_l2/bias/Read/ReadVariableOp%decoder_l3/kernel/Read/ReadVariableOp#decoder_l3/bias/Read/ReadVariableOp%decoder_l4/kernel/Read/ReadVariableOp#decoder_l4/bias/Read/ReadVariableOp'output_layer/kernel/Read/ReadVariableOp%output_layer/bias/Read/ReadVariableOptotal/Read/ReadVariableOpcount/Read/ReadVariableOp,Adam/encoder_l2/kernel/m/Read/ReadVariableOp*Adam/encoder_l2/bias/m/Read/ReadVariableOp,Adam/encoder_l3/kernel/m/Read/ReadVariableOp*Adam/encoder_l3/bias/m/Read/ReadVariableOp,Adam/encoder_l4/kernel/m/Read/ReadVariableOp*Adam/encoder_l4/bias/m/Read/ReadVariableOp(Adam/z_mean/kernel/m/Read/ReadVariableOp&Adam/z_mean/bias/m/Read/ReadVariableOp+Adam/z_log_var/kernel/m/Read/ReadVariableOp)Adam/z_log_var/bias/m/Read/ReadVariableOp,Adam/decoder_l2/kernel/m/Read/ReadVariableOp*Adam/decoder_l2/bias/m/Read/ReadVariableOp,Adam/decoder_l3/kernel/m/Read/ReadVariableOp*Adam/decoder_l3/bias/m/Read/ReadVariableOp,Adam/decoder_l4/kernel/m/Read/ReadVariableOp*Adam/decoder_l4/bias/m/Read/ReadVariableOp.Adam/output_layer/kernel/m/Read/ReadVariableOp,Adam/output_layer/bias/m/Read/ReadVariableOp,Adam/encoder_l2/kernel/v/Read/ReadVariableOp*Adam/encoder_l2/bias/v/Read/ReadVariableOp,Adam/encoder_l3/kernel/v/Read/ReadVariableOp*Adam/encoder_l3/bias/v/Read/ReadVariableOp,Adam/encoder_l4/kernel/v/Read/ReadVariableOp*Adam/encoder_l4/bias/v/Read/ReadVariableOp(Adam/z_mean/kernel/v/Read/ReadVariableOp&Adam/z_mean/bias/v/Read/ReadVariableOp+Adam/z_log_var/kernel/v/Read/ReadVariableOp)Adam/z_log_var/bias/v/Read/ReadVariableOp,Adam/decoder_l2/kernel/v/Read/ReadVariableOp*Adam/decoder_l2/bias/v/Read/ReadVariableOp,Adam/decoder_l3/kernel/v/Read/ReadVariableOp*Adam/decoder_l3/bias/v/Read/ReadVariableOp,Adam/decoder_l4/kernel/v/Read/ReadVariableOp*Adam/decoder_l4/bias/v/Read/ReadVariableOp.Adam/output_layer/kernel/v/Read/ReadVariableOp,Adam/output_layer/bias/v/Read/ReadVariableOpConst*J
TinC
A2?	*
Tout
2*
_collective_manager_ids
 *
_output_shapes
: * 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8� *(
f#R!
__inference__traced_save_143327
�
StatefulPartitionedCall_2StatefulPartitionedCallsaver_filename	Adam/iterAdam/beta_1Adam/beta_2
Adam/decayAdam/learning_rateencoder_l2/kernelencoder_l2/biasencoder_l3/kernelencoder_l3/biasencoder_l4/kernelencoder_l4/biasz_mean/kernelz_mean/biasz_log_var/kernelz_log_var/biasdecoder_l2/kerneldecoder_l2/biasdecoder_l3/kerneldecoder_l3/biasdecoder_l4/kerneldecoder_l4/biasoutput_layer/kerneloutput_layer/biastotalcountAdam/encoder_l2/kernel/mAdam/encoder_l2/bias/mAdam/encoder_l3/kernel/mAdam/encoder_l3/bias/mAdam/encoder_l4/kernel/mAdam/encoder_l4/bias/mAdam/z_mean/kernel/mAdam/z_mean/bias/mAdam/z_log_var/kernel/mAdam/z_log_var/bias/mAdam/decoder_l2/kernel/mAdam/decoder_l2/bias/mAdam/decoder_l3/kernel/mAdam/decoder_l3/bias/mAdam/decoder_l4/kernel/mAdam/decoder_l4/bias/mAdam/output_layer/kernel/mAdam/output_layer/bias/mAdam/encoder_l2/kernel/vAdam/encoder_l2/bias/vAdam/encoder_l3/kernel/vAdam/encoder_l3/bias/vAdam/encoder_l4/kernel/vAdam/encoder_l4/bias/vAdam/z_mean/kernel/vAdam/z_mean/bias/vAdam/z_log_var/kernel/vAdam/z_log_var/bias/vAdam/decoder_l2/kernel/vAdam/decoder_l2/bias/vAdam/decoder_l3/kernel/vAdam/decoder_l3/bias/vAdam/decoder_l4/kernel/vAdam/decoder_l4/bias/vAdam/output_layer/kernel/vAdam/output_layer/bias/v*I
TinB
@2>*
Tout
2*
_collective_manager_ids
 *
_output_shapes
: * 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8� *+
f&R$
"__inference__traced_restore_143520��
�
�
(__inference_encoder_layer_call_fn_141482
encoder_input
unknown
	unknown_0
	unknown_1
	unknown_2
	unknown_3
	unknown_4
	unknown_5
	unknown_6
	unknown_7
	unknown_8
identity

identity_1

identity_2��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallencoder_inputunknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6	unknown_7	unknown_8*
Tin
2*
Tout
2*
_collective_manager_ids
 *M
_output_shapes;
9:���������:���������:���������*,
_read_only_resource_inputs

	
*-
config_proto

CPU

GPU 2J 8� *L
fGRE
C__inference_encoder_layer_call_and_return_conditional_losses_1414552
StatefulPartitionedCall�
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*'
_output_shapes
:���������2

Identity�

Identity_1Identity StatefulPartitionedCall:output:1^StatefulPartitionedCall*
T0*'
_output_shapes
:���������2

Identity_1�

Identity_2Identity StatefulPartitionedCall:output:2^StatefulPartitionedCall*
T0*'
_output_shapes
:���������2

Identity_2"
identityIdentity:output:0"!

identity_1Identity_1:output:0"!

identity_2Identity_2:output:0*N
_input_shapes=
;:���������::::::::::22
StatefulPartitionedCallStatefulPartitionedCall:V R
'
_output_shapes
:���������
'
_user_specified_nameencoder_input
�	
�
E__inference_z_log_var_layer_call_and_return_conditional_losses_141273

inputs"
matmul_readvariableop_resource#
biasadd_readvariableop_resource
identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOp�
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:*
dtype02
MatMul/ReadVariableOps
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2
MatMul�
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:*
dtype02
BiasAdd/ReadVariableOp�
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2	
BiasAdd�
IdentityIdentityBiasAdd:output:0^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*
T0*'
_output_shapes
:���������2

Identity"
identityIdentity:output:0*.
_input_shapes
:���������::20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:O K
'
_output_shapes
:���������
 
_user_specified_nameinputs
�&
�
C__inference_encoder_layer_call_and_return_conditional_losses_141327
encoder_input
encoder_l2_141178
encoder_l2_141180
encoder_l3_141205
encoder_l3_141207
encoder_l4_141232
encoder_l4_141234
z_mean_141258
z_mean_141260
z_log_var_141284
z_log_var_141286
identity

identity_1

identity_2��"encoder_l2/StatefulPartitionedCall�"encoder_l3/StatefulPartitionedCall�"encoder_l4/StatefulPartitionedCall�"sampling_3/StatefulPartitionedCall�!z_log_var/StatefulPartitionedCall�z_mean/StatefulPartitionedCall�
"encoder_l2/StatefulPartitionedCallStatefulPartitionedCallencoder_inputencoder_l2_141178encoder_l2_141180*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������
*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *O
fJRH
F__inference_encoder_l2_layer_call_and_return_conditional_losses_1411672$
"encoder_l2/StatefulPartitionedCall�
"encoder_l3/StatefulPartitionedCallStatefulPartitionedCall+encoder_l2/StatefulPartitionedCall:output:0encoder_l3_141205encoder_l3_141207*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *O
fJRH
F__inference_encoder_l3_layer_call_and_return_conditional_losses_1411942$
"encoder_l3/StatefulPartitionedCall�
"encoder_l4/StatefulPartitionedCallStatefulPartitionedCall+encoder_l3/StatefulPartitionedCall:output:0encoder_l4_141232encoder_l4_141234*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *O
fJRH
F__inference_encoder_l4_layer_call_and_return_conditional_losses_1412212$
"encoder_l4/StatefulPartitionedCall�
z_mean/StatefulPartitionedCallStatefulPartitionedCall+encoder_l4/StatefulPartitionedCall:output:0z_mean_141258z_mean_141260*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *K
fFRD
B__inference_z_mean_layer_call_and_return_conditional_losses_1412472 
z_mean/StatefulPartitionedCall�
!z_log_var/StatefulPartitionedCallStatefulPartitionedCall+encoder_l4/StatefulPartitionedCall:output:0z_log_var_141284z_log_var_141286*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *N
fIRG
E__inference_z_log_var_layer_call_and_return_conditional_losses_1412732#
!z_log_var/StatefulPartitionedCall�
"sampling_3/StatefulPartitionedCallStatefulPartitionedCall'z_mean/StatefulPartitionedCall:output:0*z_log_var/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8� *O
fJRH
F__inference_sampling_3_layer_call_and_return_conditional_losses_1413152$
"sampling_3/StatefulPartitionedCall�
IdentityIdentity'z_mean/StatefulPartitionedCall:output:0#^encoder_l2/StatefulPartitionedCall#^encoder_l3/StatefulPartitionedCall#^encoder_l4/StatefulPartitionedCall#^sampling_3/StatefulPartitionedCall"^z_log_var/StatefulPartitionedCall^z_mean/StatefulPartitionedCall*
T0*'
_output_shapes
:���������2

Identity�

Identity_1Identity*z_log_var/StatefulPartitionedCall:output:0#^encoder_l2/StatefulPartitionedCall#^encoder_l3/StatefulPartitionedCall#^encoder_l4/StatefulPartitionedCall#^sampling_3/StatefulPartitionedCall"^z_log_var/StatefulPartitionedCall^z_mean/StatefulPartitionedCall*
T0*'
_output_shapes
:���������2

Identity_1�

Identity_2Identity+sampling_3/StatefulPartitionedCall:output:0#^encoder_l2/StatefulPartitionedCall#^encoder_l3/StatefulPartitionedCall#^encoder_l4/StatefulPartitionedCall#^sampling_3/StatefulPartitionedCall"^z_log_var/StatefulPartitionedCall^z_mean/StatefulPartitionedCall*
T0*'
_output_shapes
:���������2

Identity_2"
identityIdentity:output:0"!

identity_1Identity_1:output:0"!

identity_2Identity_2:output:0*N
_input_shapes=
;:���������::::::::::2H
"encoder_l2/StatefulPartitionedCall"encoder_l2/StatefulPartitionedCall2H
"encoder_l3/StatefulPartitionedCall"encoder_l3/StatefulPartitionedCall2H
"encoder_l4/StatefulPartitionedCall"encoder_l4/StatefulPartitionedCall2H
"sampling_3/StatefulPartitionedCall"sampling_3/StatefulPartitionedCall2F
!z_log_var/StatefulPartitionedCall!z_log_var/StatefulPartitionedCall2@
z_mean/StatefulPartitionedCallz_mean/StatefulPartitionedCall:V R
'
_output_shapes
:���������
'
_user_specified_nameencoder_input
�
�
(__inference_decoder_layer_call_fn_141710
latent_variable
unknown
	unknown_0
	unknown_1
	unknown_2
	unknown_3
	unknown_4
	unknown_5
	unknown_6
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCalllatent_variableunknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6*
Tin
2	*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������**
_read_only_resource_inputs

*-
config_proto

CPU

GPU 2J 8� *L
fGRE
C__inference_decoder_layer_call_and_return_conditional_losses_1416912
StatefulPartitionedCall�
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*'
_output_shapes
:���������2

Identity"
identityIdentity:output:0*F
_input_shapes5
3:���������::::::::22
StatefulPartitionedCallStatefulPartitionedCall:X T
'
_output_shapes
:���������
)
_user_specified_namelatent_variable
�
�
+__inference_encoder_l3_layer_call_fn_142951

inputs
unknown
	unknown_0
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *O
fJRH
F__inference_encoder_l3_layer_call_and_return_conditional_losses_1411942
StatefulPartitionedCall�
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*'
_output_shapes
:���������2

Identity"
identityIdentity:output:0*.
_input_shapes
:���������
::22
StatefulPartitionedCallStatefulPartitionedCall:O K
'
_output_shapes
:���������

 
_user_specified_nameinputs
�
�
C__inference_decoder_layer_call_and_return_conditional_losses_141595
latent_variable
decoder_l2_141508
decoder_l2_141510
decoder_l3_141535
decoder_l3_141537
decoder_l4_141562
decoder_l4_141564
output_layer_141589
output_layer_141591
identity��"decoder_l2/StatefulPartitionedCall�"decoder_l3/StatefulPartitionedCall�"decoder_l4/StatefulPartitionedCall�$output_layer/StatefulPartitionedCall�
"decoder_l2/StatefulPartitionedCallStatefulPartitionedCalllatent_variabledecoder_l2_141508decoder_l2_141510*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *O
fJRH
F__inference_decoder_l2_layer_call_and_return_conditional_losses_1414972$
"decoder_l2/StatefulPartitionedCall�
"decoder_l3/StatefulPartitionedCallStatefulPartitionedCall+decoder_l2/StatefulPartitionedCall:output:0decoder_l3_141535decoder_l3_141537*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *O
fJRH
F__inference_decoder_l3_layer_call_and_return_conditional_losses_1415242$
"decoder_l3/StatefulPartitionedCall�
"decoder_l4/StatefulPartitionedCallStatefulPartitionedCall+decoder_l3/StatefulPartitionedCall:output:0decoder_l4_141562decoder_l4_141564*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������
*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *O
fJRH
F__inference_decoder_l4_layer_call_and_return_conditional_losses_1415512$
"decoder_l4/StatefulPartitionedCall�
$output_layer/StatefulPartitionedCallStatefulPartitionedCall+decoder_l4/StatefulPartitionedCall:output:0output_layer_141589output_layer_141591*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *Q
fLRJ
H__inference_output_layer_layer_call_and_return_conditional_losses_1415782&
$output_layer/StatefulPartitionedCall�
IdentityIdentity-output_layer/StatefulPartitionedCall:output:0#^decoder_l2/StatefulPartitionedCall#^decoder_l3/StatefulPartitionedCall#^decoder_l4/StatefulPartitionedCall%^output_layer/StatefulPartitionedCall*
T0*'
_output_shapes
:���������2

Identity"
identityIdentity:output:0*F
_input_shapes5
3:���������::::::::2H
"decoder_l2/StatefulPartitionedCall"decoder_l2/StatefulPartitionedCall2H
"decoder_l3/StatefulPartitionedCall"decoder_l3/StatefulPartitionedCall2H
"decoder_l4/StatefulPartitionedCall"decoder_l4/StatefulPartitionedCall2L
$output_layer/StatefulPartitionedCall$output_layer/StatefulPartitionedCall:X T
'
_output_shapes
:���������
)
_user_specified_namelatent_variable
�
�
(__inference_encoder_layer_call_fn_142776

inputs
unknown
	unknown_0
	unknown_1
	unknown_2
	unknown_3
	unknown_4
	unknown_5
	unknown_6
	unknown_7
	unknown_8
identity

identity_1

identity_2��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6	unknown_7	unknown_8*
Tin
2*
Tout
2*
_collective_manager_ids
 *M
_output_shapes;
9:���������:���������:���������*,
_read_only_resource_inputs

	
*-
config_proto

CPU

GPU 2J 8� *L
fGRE
C__inference_encoder_layer_call_and_return_conditional_losses_1413942
StatefulPartitionedCall�
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*'
_output_shapes
:���������2

Identity�

Identity_1Identity StatefulPartitionedCall:output:1^StatefulPartitionedCall*
T0*'
_output_shapes
:���������2

Identity_1�

Identity_2Identity StatefulPartitionedCall:output:2^StatefulPartitionedCall*
T0*'
_output_shapes
:���������2

Identity_2"
identityIdentity:output:0"!

identity_1Identity_1:output:0"!

identity_2Identity_2:output:0*N
_input_shapes=
;:���������::::::::::22
StatefulPartitionedCallStatefulPartitionedCall:O K
'
_output_shapes
:���������
 
_user_specified_nameinputs
�	
�
F__inference_encoder_l4_layer_call_and_return_conditional_losses_142962

inputs"
matmul_readvariableop_resource#
biasadd_readvariableop_resource
identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOp�
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:*
dtype02
MatMul/ReadVariableOps
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2
MatMul�
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:*
dtype02
BiasAdd/ReadVariableOp�
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2	
BiasAddX
ReluReluBiasAdd:output:0*
T0*'
_output_shapes
:���������2
Relu�
IdentityIdentityRelu:activations:0^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*
T0*'
_output_shapes
:���������2

Identity"
identityIdentity:output:0*.
_input_shapes
:���������::20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:O K
'
_output_shapes
:���������
 
_user_specified_nameinputs
�	
�
E__inference_z_log_var_layer_call_and_return_conditional_losses_143000

inputs"
matmul_readvariableop_resource#
biasadd_readvariableop_resource
identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOp�
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:*
dtype02
MatMul/ReadVariableOps
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2
MatMul�
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:*
dtype02
BiasAdd/ReadVariableOp�
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2	
BiasAdd�
IdentityIdentityBiasAdd:output:0^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*
T0*'
_output_shapes
:���������2

Identity"
identityIdentity:output:0*.
_input_shapes
:���������::20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:O K
'
_output_shapes
:���������
 
_user_specified_nameinputs
�
�
+__inference_decoder_l3_layer_call_fn_143081

inputs
unknown
	unknown_0
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *O
fJRH
F__inference_decoder_l3_layer_call_and_return_conditional_losses_1415242
StatefulPartitionedCall�
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*'
_output_shapes
:���������2

Identity"
identityIdentity:output:0*.
_input_shapes
:���������::22
StatefulPartitionedCallStatefulPartitionedCall:O K
'
_output_shapes
:���������
 
_user_specified_nameinputs
�R
�
C__inference_encoder_layer_call_and_return_conditional_losses_142687

inputs-
)encoder_l2_matmul_readvariableop_resource.
*encoder_l2_biasadd_readvariableop_resource-
)encoder_l3_matmul_readvariableop_resource.
*encoder_l3_biasadd_readvariableop_resource-
)encoder_l4_matmul_readvariableop_resource.
*encoder_l4_biasadd_readvariableop_resource)
%z_mean_matmul_readvariableop_resource*
&z_mean_biasadd_readvariableop_resource,
(z_log_var_matmul_readvariableop_resource-
)z_log_var_biasadd_readvariableop_resource
identity

identity_1

identity_2��!encoder_l2/BiasAdd/ReadVariableOp� encoder_l2/MatMul/ReadVariableOp�!encoder_l3/BiasAdd/ReadVariableOp� encoder_l3/MatMul/ReadVariableOp�!encoder_l4/BiasAdd/ReadVariableOp� encoder_l4/MatMul/ReadVariableOp� z_log_var/BiasAdd/ReadVariableOp�z_log_var/MatMul/ReadVariableOp�z_mean/BiasAdd/ReadVariableOp�z_mean/MatMul/ReadVariableOp�
 encoder_l2/MatMul/ReadVariableOpReadVariableOp)encoder_l2_matmul_readvariableop_resource*
_output_shapes

:
*
dtype02"
 encoder_l2/MatMul/ReadVariableOp�
encoder_l2/MatMulMatMulinputs(encoder_l2/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������
2
encoder_l2/MatMul�
!encoder_l2/BiasAdd/ReadVariableOpReadVariableOp*encoder_l2_biasadd_readvariableop_resource*
_output_shapes
:
*
dtype02#
!encoder_l2/BiasAdd/ReadVariableOp�
encoder_l2/BiasAddBiasAddencoder_l2/MatMul:product:0)encoder_l2/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������
2
encoder_l2/BiasAddy
encoder_l2/ReluReluencoder_l2/BiasAdd:output:0*
T0*'
_output_shapes
:���������
2
encoder_l2/Relu�
 encoder_l3/MatMul/ReadVariableOpReadVariableOp)encoder_l3_matmul_readvariableop_resource*
_output_shapes

:
*
dtype02"
 encoder_l3/MatMul/ReadVariableOp�
encoder_l3/MatMulMatMulencoder_l2/Relu:activations:0(encoder_l3/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2
encoder_l3/MatMul�
!encoder_l3/BiasAdd/ReadVariableOpReadVariableOp*encoder_l3_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02#
!encoder_l3/BiasAdd/ReadVariableOp�
encoder_l3/BiasAddBiasAddencoder_l3/MatMul:product:0)encoder_l3/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2
encoder_l3/BiasAddy
encoder_l3/ReluReluencoder_l3/BiasAdd:output:0*
T0*'
_output_shapes
:���������2
encoder_l3/Relu�
 encoder_l4/MatMul/ReadVariableOpReadVariableOp)encoder_l4_matmul_readvariableop_resource*
_output_shapes

:*
dtype02"
 encoder_l4/MatMul/ReadVariableOp�
encoder_l4/MatMulMatMulencoder_l3/Relu:activations:0(encoder_l4/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2
encoder_l4/MatMul�
!encoder_l4/BiasAdd/ReadVariableOpReadVariableOp*encoder_l4_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02#
!encoder_l4/BiasAdd/ReadVariableOp�
encoder_l4/BiasAddBiasAddencoder_l4/MatMul:product:0)encoder_l4/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2
encoder_l4/BiasAddy
encoder_l4/ReluReluencoder_l4/BiasAdd:output:0*
T0*'
_output_shapes
:���������2
encoder_l4/Relu�
z_mean/MatMul/ReadVariableOpReadVariableOp%z_mean_matmul_readvariableop_resource*
_output_shapes

:*
dtype02
z_mean/MatMul/ReadVariableOp�
z_mean/MatMulMatMulencoder_l4/Relu:activations:0$z_mean/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2
z_mean/MatMul�
z_mean/BiasAdd/ReadVariableOpReadVariableOp&z_mean_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02
z_mean/BiasAdd/ReadVariableOp�
z_mean/BiasAddBiasAddz_mean/MatMul:product:0%z_mean/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2
z_mean/BiasAdd�
z_log_var/MatMul/ReadVariableOpReadVariableOp(z_log_var_matmul_readvariableop_resource*
_output_shapes

:*
dtype02!
z_log_var/MatMul/ReadVariableOp�
z_log_var/MatMulMatMulencoder_l4/Relu:activations:0'z_log_var/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2
z_log_var/MatMul�
 z_log_var/BiasAdd/ReadVariableOpReadVariableOp)z_log_var_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02"
 z_log_var/BiasAdd/ReadVariableOp�
z_log_var/BiasAddBiasAddz_log_var/MatMul:product:0(z_log_var/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2
z_log_var/BiasAddk
sampling_3/ShapeShapez_mean/BiasAdd:output:0*
T0*
_output_shapes
:2
sampling_3/Shape�
sampling_3/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB: 2 
sampling_3/strided_slice/stack�
 sampling_3/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2"
 sampling_3/strided_slice/stack_1�
 sampling_3/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2"
 sampling_3/strided_slice/stack_2�
sampling_3/strided_sliceStridedSlicesampling_3/Shape:output:0'sampling_3/strided_slice/stack:output:0)sampling_3/strided_slice/stack_1:output:0)sampling_3/strided_slice/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2
sampling_3/strided_sliceo
sampling_3/Shape_1Shapez_mean/BiasAdd:output:0*
T0*
_output_shapes
:2
sampling_3/Shape_1�
 sampling_3/strided_slice_1/stackConst*
_output_shapes
:*
dtype0*
valueB:2"
 sampling_3/strided_slice_1/stack�
"sampling_3/strided_slice_1/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2$
"sampling_3/strided_slice_1/stack_1�
"sampling_3/strided_slice_1/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2$
"sampling_3/strided_slice_1/stack_2�
sampling_3/strided_slice_1StridedSlicesampling_3/Shape_1:output:0)sampling_3/strided_slice_1/stack:output:0+sampling_3/strided_slice_1/stack_1:output:0+sampling_3/strided_slice_1/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2
sampling_3/strided_slice_1�
sampling_3/random_normal/shapePack!sampling_3/strided_slice:output:0#sampling_3/strided_slice_1:output:0*
N*
T0*
_output_shapes
:2 
sampling_3/random_normal/shape�
sampling_3/random_normal/meanConst*
_output_shapes
: *
dtype0*
valueB
 *    2
sampling_3/random_normal/mean�
sampling_3/random_normal/stddevConst*
_output_shapes
: *
dtype0*
valueB
 *  �?2!
sampling_3/random_normal/stddev�
-sampling_3/random_normal/RandomStandardNormalRandomStandardNormal'sampling_3/random_normal/shape:output:0*
T0*0
_output_shapes
:������������������*
dtype0*
seed���)*
seed2Đ�2/
-sampling_3/random_normal/RandomStandardNormal�
sampling_3/random_normal/mulMul6sampling_3/random_normal/RandomStandardNormal:output:0(sampling_3/random_normal/stddev:output:0*
T0*0
_output_shapes
:������������������2
sampling_3/random_normal/mul�
sampling_3/random_normalAdd sampling_3/random_normal/mul:z:0&sampling_3/random_normal/mean:output:0*
T0*0
_output_shapes
:������������������2
sampling_3/random_normali
sampling_3/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *   ?2
sampling_3/mul/x�
sampling_3/mulMulsampling_3/mul/x:output:0z_log_var/BiasAdd:output:0*
T0*'
_output_shapes
:���������2
sampling_3/mulm
sampling_3/ExpExpsampling_3/mul:z:0*
T0*'
_output_shapes
:���������2
sampling_3/Exp�
sampling_3/mul_1Mulsampling_3/Exp:y:0sampling_3/random_normal:z:0*
T0*'
_output_shapes
:���������2
sampling_3/mul_1�
sampling_3/addAddV2z_mean/BiasAdd:output:0sampling_3/mul_1:z:0*
T0*'
_output_shapes
:���������2
sampling_3/add�
IdentityIdentityz_mean/BiasAdd:output:0"^encoder_l2/BiasAdd/ReadVariableOp!^encoder_l2/MatMul/ReadVariableOp"^encoder_l3/BiasAdd/ReadVariableOp!^encoder_l3/MatMul/ReadVariableOp"^encoder_l4/BiasAdd/ReadVariableOp!^encoder_l4/MatMul/ReadVariableOp!^z_log_var/BiasAdd/ReadVariableOp ^z_log_var/MatMul/ReadVariableOp^z_mean/BiasAdd/ReadVariableOp^z_mean/MatMul/ReadVariableOp*
T0*'
_output_shapes
:���������2

Identity�

Identity_1Identityz_log_var/BiasAdd:output:0"^encoder_l2/BiasAdd/ReadVariableOp!^encoder_l2/MatMul/ReadVariableOp"^encoder_l3/BiasAdd/ReadVariableOp!^encoder_l3/MatMul/ReadVariableOp"^encoder_l4/BiasAdd/ReadVariableOp!^encoder_l4/MatMul/ReadVariableOp!^z_log_var/BiasAdd/ReadVariableOp ^z_log_var/MatMul/ReadVariableOp^z_mean/BiasAdd/ReadVariableOp^z_mean/MatMul/ReadVariableOp*
T0*'
_output_shapes
:���������2

Identity_1�

Identity_2Identitysampling_3/add:z:0"^encoder_l2/BiasAdd/ReadVariableOp!^encoder_l2/MatMul/ReadVariableOp"^encoder_l3/BiasAdd/ReadVariableOp!^encoder_l3/MatMul/ReadVariableOp"^encoder_l4/BiasAdd/ReadVariableOp!^encoder_l4/MatMul/ReadVariableOp!^z_log_var/BiasAdd/ReadVariableOp ^z_log_var/MatMul/ReadVariableOp^z_mean/BiasAdd/ReadVariableOp^z_mean/MatMul/ReadVariableOp*
T0*'
_output_shapes
:���������2

Identity_2"
identityIdentity:output:0"!

identity_1Identity_1:output:0"!

identity_2Identity_2:output:0*N
_input_shapes=
;:���������::::::::::2F
!encoder_l2/BiasAdd/ReadVariableOp!encoder_l2/BiasAdd/ReadVariableOp2D
 encoder_l2/MatMul/ReadVariableOp encoder_l2/MatMul/ReadVariableOp2F
!encoder_l3/BiasAdd/ReadVariableOp!encoder_l3/BiasAdd/ReadVariableOp2D
 encoder_l3/MatMul/ReadVariableOp encoder_l3/MatMul/ReadVariableOp2F
!encoder_l4/BiasAdd/ReadVariableOp!encoder_l4/BiasAdd/ReadVariableOp2D
 encoder_l4/MatMul/ReadVariableOp encoder_l4/MatMul/ReadVariableOp2D
 z_log_var/BiasAdd/ReadVariableOp z_log_var/BiasAdd/ReadVariableOp2B
z_log_var/MatMul/ReadVariableOpz_log_var/MatMul/ReadVariableOp2>
z_mean/BiasAdd/ReadVariableOpz_mean/BiasAdd/ReadVariableOp2<
z_mean/MatMul/ReadVariableOpz_mean/MatMul/ReadVariableOp:O K
'
_output_shapes
:���������
 
_user_specified_nameinputs
�
t
+__inference_sampling_3_layer_call_fn_143041
inputs_0
inputs_1
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputs_0inputs_1*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8� *O
fJRH
F__inference_sampling_3_layer_call_and_return_conditional_losses_1413152
StatefulPartitionedCall�
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*'
_output_shapes
:���������2

Identity"
identityIdentity:output:0*9
_input_shapes(
&:���������:���������22
StatefulPartitionedCallStatefulPartitionedCall:Q M
'
_output_shapes
:���������
"
_user_specified_name
inputs/0:QM
'
_output_shapes
:���������
"
_user_specified_name
inputs/1
�
�
(__inference_decoder_layer_call_fn_142890

inputs
unknown
	unknown_0
	unknown_1
	unknown_2
	unknown_3
	unknown_4
	unknown_5
	unknown_6
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6*
Tin
2	*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������**
_read_only_resource_inputs

*-
config_proto

CPU

GPU 2J 8� *L
fGRE
C__inference_decoder_layer_call_and_return_conditional_losses_1416462
StatefulPartitionedCall�
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*'
_output_shapes
:���������2

Identity"
identityIdentity:output:0*F
_input_shapes5
3:���������::::::::22
StatefulPartitionedCallStatefulPartitionedCall:O K
'
_output_shapes
:���������
 
_user_specified_nameinputs
�
�
B__inference_cvae_3_layer_call_and_return_conditional_losses_141964

inputs
encoder_141917
encoder_141919
encoder_141921
encoder_141923
encoder_141925
encoder_141927
encoder_141929
encoder_141931
encoder_141933
encoder_141935
decoder_141946
decoder_141948
decoder_141950
decoder_141952
decoder_141954
decoder_141956
decoder_141958
decoder_141960
identity��decoder/StatefulPartitionedCall�encoder/StatefulPartitionedCall�
encoder/StatefulPartitionedCallStatefulPartitionedCallinputsencoder_141917encoder_141919encoder_141921encoder_141923encoder_141925encoder_141927encoder_141929encoder_141931encoder_141933encoder_141935*
Tin
2*
Tout
2*
_collective_manager_ids
 *M
_output_shapes;
9:���������:���������:���������*,
_read_only_resource_inputs

	
*-
config_proto

CPU

GPU 2J 8� *L
fGRE
C__inference_encoder_layer_call_and_return_conditional_losses_1414552!
encoder/StatefulPartitionedCall{
strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB"    ����2
strided_slice/stack
strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB"        2
strided_slice/stack_1
strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      2
strided_slice/stack_2�
strided_sliceStridedSliceinputsstrided_slice/stack:output:0strided_slice/stack_1:output:0strided_slice/stack_2:output:0*
Index0*
T0*'
_output_shapes
:���������*

begin_mask*
end_mask2
strided_slicet
concatenate/concat/axisConst*
_output_shapes
: *
dtype0*
value	B :2
concatenate/concat/axis�
concatenate/concatConcatV2(encoder/StatefulPartitionedCall:output:2strided_slice:output:0 concatenate/concat/axis:output:0*
N*
T0*'
_output_shapes
:���������2
concatenate/concat�
decoder/StatefulPartitionedCallStatefulPartitionedCallconcatenate/concat:output:0decoder_141946decoder_141948decoder_141950decoder_141952decoder_141954decoder_141956decoder_141958decoder_141960*
Tin
2	*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������**
_read_only_resource_inputs

*-
config_proto

CPU

GPU 2J 8� *L
fGRE
C__inference_decoder_layer_call_and_return_conditional_losses_1416912!
decoder/StatefulPartitionedCall�
IdentityIdentity(decoder/StatefulPartitionedCall:output:0 ^decoder/StatefulPartitionedCall ^encoder/StatefulPartitionedCall*
T0*'
_output_shapes
:���������2

Identity"
identityIdentity:output:0*n
_input_shapes]
[:���������::::::::::::::::::2B
decoder/StatefulPartitionedCalldecoder/StatefulPartitionedCall2B
encoder/StatefulPartitionedCallencoder/StatefulPartitionedCall:O K
'
_output_shapes
:���������
 
_user_specified_nameinputs
�
�
C__inference_decoder_layer_call_and_return_conditional_losses_141619
latent_variable
decoder_l2_141598
decoder_l2_141600
decoder_l3_141603
decoder_l3_141605
decoder_l4_141608
decoder_l4_141610
output_layer_141613
output_layer_141615
identity��"decoder_l2/StatefulPartitionedCall�"decoder_l3/StatefulPartitionedCall�"decoder_l4/StatefulPartitionedCall�$output_layer/StatefulPartitionedCall�
"decoder_l2/StatefulPartitionedCallStatefulPartitionedCalllatent_variabledecoder_l2_141598decoder_l2_141600*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *O
fJRH
F__inference_decoder_l2_layer_call_and_return_conditional_losses_1414972$
"decoder_l2/StatefulPartitionedCall�
"decoder_l3/StatefulPartitionedCallStatefulPartitionedCall+decoder_l2/StatefulPartitionedCall:output:0decoder_l3_141603decoder_l3_141605*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *O
fJRH
F__inference_decoder_l3_layer_call_and_return_conditional_losses_1415242$
"decoder_l3/StatefulPartitionedCall�
"decoder_l4/StatefulPartitionedCallStatefulPartitionedCall+decoder_l3/StatefulPartitionedCall:output:0decoder_l4_141608decoder_l4_141610*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������
*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *O
fJRH
F__inference_decoder_l4_layer_call_and_return_conditional_losses_1415512$
"decoder_l4/StatefulPartitionedCall�
$output_layer/StatefulPartitionedCallStatefulPartitionedCall+decoder_l4/StatefulPartitionedCall:output:0output_layer_141613output_layer_141615*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *Q
fLRJ
H__inference_output_layer_layer_call_and_return_conditional_losses_1415782&
$output_layer/StatefulPartitionedCall�
IdentityIdentity-output_layer/StatefulPartitionedCall:output:0#^decoder_l2/StatefulPartitionedCall#^decoder_l3/StatefulPartitionedCall#^decoder_l4/StatefulPartitionedCall%^output_layer/StatefulPartitionedCall*
T0*'
_output_shapes
:���������2

Identity"
identityIdentity:output:0*F
_input_shapes5
3:���������::::::::2H
"decoder_l2/StatefulPartitionedCall"decoder_l2/StatefulPartitionedCall2H
"decoder_l3/StatefulPartitionedCall"decoder_l3/StatefulPartitionedCall2H
"decoder_l4/StatefulPartitionedCall"decoder_l4/StatefulPartitionedCall2L
$output_layer/StatefulPartitionedCall$output_layer/StatefulPartitionedCall:X T
'
_output_shapes
:���������
)
_user_specified_namelatent_variable
�
�
+__inference_decoder_l4_layer_call_fn_143101

inputs
unknown
	unknown_0
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������
*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *O
fJRH
F__inference_decoder_l4_layer_call_and_return_conditional_losses_1415512
StatefulPartitionedCall�
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*'
_output_shapes
:���������
2

Identity"
identityIdentity:output:0*.
_input_shapes
:���������::22
StatefulPartitionedCallStatefulPartitionedCall:O K
'
_output_shapes
:���������
 
_user_specified_nameinputs
�
�
+__inference_decoder_l2_layer_call_fn_143061

inputs
unknown
	unknown_0
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *O
fJRH
F__inference_decoder_l2_layer_call_and_return_conditional_losses_1414972
StatefulPartitionedCall�
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*'
_output_shapes
:���������2

Identity"
identityIdentity:output:0*.
_input_shapes
:���������::22
StatefulPartitionedCallStatefulPartitionedCall:O K
'
_output_shapes
:���������
 
_user_specified_nameinputs
�	
�
H__inference_output_layer_layer_call_and_return_conditional_losses_141578

inputs"
matmul_readvariableop_resource#
biasadd_readvariableop_resource
identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOp�
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:
*
dtype02
MatMul/ReadVariableOps
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2
MatMul�
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:*
dtype02
BiasAdd/ReadVariableOp�
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2	
BiasAddX
ReluReluBiasAdd:output:0*
T0*'
_output_shapes
:���������2
Relu�
IdentityIdentityRelu:activations:0^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*
T0*'
_output_shapes
:���������2

Identity"
identityIdentity:output:0*.
_input_shapes
:���������
::20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:O K
'
_output_shapes
:���������

 
_user_specified_nameinputs
�	
�
F__inference_encoder_l2_layer_call_and_return_conditional_losses_141167

inputs"
matmul_readvariableop_resource#
biasadd_readvariableop_resource
identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOp�
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:
*
dtype02
MatMul/ReadVariableOps
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������
2
MatMul�
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:
*
dtype02
BiasAdd/ReadVariableOp�
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������
2	
BiasAddX
ReluReluBiasAdd:output:0*
T0*'
_output_shapes
:���������
2
Relu�
IdentityIdentityRelu:activations:0^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*
T0*'
_output_shapes
:���������
2

Identity"
identityIdentity:output:0*.
_input_shapes
:���������::20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:O K
'
_output_shapes
:���������
 
_user_specified_nameinputs
�	
�
F__inference_decoder_l4_layer_call_and_return_conditional_losses_143092

inputs"
matmul_readvariableop_resource#
biasadd_readvariableop_resource
identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOp�
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:
*
dtype02
MatMul/ReadVariableOps
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������
2
MatMul�
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:
*
dtype02
BiasAdd/ReadVariableOp�
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������
2	
BiasAddX
ReluReluBiasAdd:output:0*
T0*'
_output_shapes
:���������
2
Relu�
IdentityIdentityRelu:activations:0^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*
T0*'
_output_shapes
:���������
2

Identity"
identityIdentity:output:0*.
_input_shapes
:���������::20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:O K
'
_output_shapes
:���������
 
_user_specified_nameinputs
�&
�
C__inference_encoder_layer_call_and_return_conditional_losses_141394

inputs
encoder_l2_141365
encoder_l2_141367
encoder_l3_141370
encoder_l3_141372
encoder_l4_141375
encoder_l4_141377
z_mean_141380
z_mean_141382
z_log_var_141385
z_log_var_141387
identity

identity_1

identity_2��"encoder_l2/StatefulPartitionedCall�"encoder_l3/StatefulPartitionedCall�"encoder_l4/StatefulPartitionedCall�"sampling_3/StatefulPartitionedCall�!z_log_var/StatefulPartitionedCall�z_mean/StatefulPartitionedCall�
"encoder_l2/StatefulPartitionedCallStatefulPartitionedCallinputsencoder_l2_141365encoder_l2_141367*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������
*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *O
fJRH
F__inference_encoder_l2_layer_call_and_return_conditional_losses_1411672$
"encoder_l2/StatefulPartitionedCall�
"encoder_l3/StatefulPartitionedCallStatefulPartitionedCall+encoder_l2/StatefulPartitionedCall:output:0encoder_l3_141370encoder_l3_141372*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *O
fJRH
F__inference_encoder_l3_layer_call_and_return_conditional_losses_1411942$
"encoder_l3/StatefulPartitionedCall�
"encoder_l4/StatefulPartitionedCallStatefulPartitionedCall+encoder_l3/StatefulPartitionedCall:output:0encoder_l4_141375encoder_l4_141377*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *O
fJRH
F__inference_encoder_l4_layer_call_and_return_conditional_losses_1412212$
"encoder_l4/StatefulPartitionedCall�
z_mean/StatefulPartitionedCallStatefulPartitionedCall+encoder_l4/StatefulPartitionedCall:output:0z_mean_141380z_mean_141382*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *K
fFRD
B__inference_z_mean_layer_call_and_return_conditional_losses_1412472 
z_mean/StatefulPartitionedCall�
!z_log_var/StatefulPartitionedCallStatefulPartitionedCall+encoder_l4/StatefulPartitionedCall:output:0z_log_var_141385z_log_var_141387*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *N
fIRG
E__inference_z_log_var_layer_call_and_return_conditional_losses_1412732#
!z_log_var/StatefulPartitionedCall�
"sampling_3/StatefulPartitionedCallStatefulPartitionedCall'z_mean/StatefulPartitionedCall:output:0*z_log_var/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8� *O
fJRH
F__inference_sampling_3_layer_call_and_return_conditional_losses_1413152$
"sampling_3/StatefulPartitionedCall�
IdentityIdentity'z_mean/StatefulPartitionedCall:output:0#^encoder_l2/StatefulPartitionedCall#^encoder_l3/StatefulPartitionedCall#^encoder_l4/StatefulPartitionedCall#^sampling_3/StatefulPartitionedCall"^z_log_var/StatefulPartitionedCall^z_mean/StatefulPartitionedCall*
T0*'
_output_shapes
:���������2

Identity�

Identity_1Identity*z_log_var/StatefulPartitionedCall:output:0#^encoder_l2/StatefulPartitionedCall#^encoder_l3/StatefulPartitionedCall#^encoder_l4/StatefulPartitionedCall#^sampling_3/StatefulPartitionedCall"^z_log_var/StatefulPartitionedCall^z_mean/StatefulPartitionedCall*
T0*'
_output_shapes
:���������2

Identity_1�

Identity_2Identity+sampling_3/StatefulPartitionedCall:output:0#^encoder_l2/StatefulPartitionedCall#^encoder_l3/StatefulPartitionedCall#^encoder_l4/StatefulPartitionedCall#^sampling_3/StatefulPartitionedCall"^z_log_var/StatefulPartitionedCall^z_mean/StatefulPartitionedCall*
T0*'
_output_shapes
:���������2

Identity_2"
identityIdentity:output:0"!

identity_1Identity_1:output:0"!

identity_2Identity_2:output:0*N
_input_shapes=
;:���������::::::::::2H
"encoder_l2/StatefulPartitionedCall"encoder_l2/StatefulPartitionedCall2H
"encoder_l3/StatefulPartitionedCall"encoder_l3/StatefulPartitionedCall2H
"encoder_l4/StatefulPartitionedCall"encoder_l4/StatefulPartitionedCall2H
"sampling_3/StatefulPartitionedCall"sampling_3/StatefulPartitionedCall2F
!z_log_var/StatefulPartitionedCall!z_log_var/StatefulPartitionedCall2@
z_mean/StatefulPartitionedCallz_mean/StatefulPartitionedCall:O K
'
_output_shapes
:���������
 
_user_specified_nameinputs
�	
�
H__inference_output_layer_layer_call_and_return_conditional_losses_143112

inputs"
matmul_readvariableop_resource#
biasadd_readvariableop_resource
identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOp�
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:
*
dtype02
MatMul/ReadVariableOps
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2
MatMul�
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:*
dtype02
BiasAdd/ReadVariableOp�
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2	
BiasAddX
ReluReluBiasAdd:output:0*
T0*'
_output_shapes
:���������2
Relu�
IdentityIdentityRelu:activations:0^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*
T0*'
_output_shapes
:���������2

Identity"
identityIdentity:output:0*.
_input_shapes
:���������
::20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:O K
'
_output_shapes
:���������

 
_user_specified_nameinputs
�	
�
F__inference_decoder_l3_layer_call_and_return_conditional_losses_143072

inputs"
matmul_readvariableop_resource#
biasadd_readvariableop_resource
identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOp�
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:*
dtype02
MatMul/ReadVariableOps
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2
MatMul�
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:*
dtype02
BiasAdd/ReadVariableOp�
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2	
BiasAddX
ReluReluBiasAdd:output:0*
T0*'
_output_shapes
:���������2
Relu�
IdentityIdentityRelu:activations:0^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*
T0*'
_output_shapes
:���������2

Identity"
identityIdentity:output:0*.
_input_shapes
:���������::20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:O K
'
_output_shapes
:���������
 
_user_specified_nameinputs
�
�
+__inference_encoder_l4_layer_call_fn_142971

inputs
unknown
	unknown_0
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *O
fJRH
F__inference_encoder_l4_layer_call_and_return_conditional_losses_1412212
StatefulPartitionedCall�
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*'
_output_shapes
:���������2

Identity"
identityIdentity:output:0*.
_input_shapes
:���������::22
StatefulPartitionedCallStatefulPartitionedCall:O K
'
_output_shapes
:���������
 
_user_specified_nameinputs
�	
�
F__inference_decoder_l2_layer_call_and_return_conditional_losses_141497

inputs"
matmul_readvariableop_resource#
biasadd_readvariableop_resource
identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOp�
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:*
dtype02
MatMul/ReadVariableOps
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2
MatMul�
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:*
dtype02
BiasAdd/ReadVariableOp�
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2	
BiasAddX
ReluReluBiasAdd:output:0*
T0*'
_output_shapes
:���������2
Relu�
IdentityIdentityRelu:activations:0^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*
T0*'
_output_shapes
:���������2

Identity"
identityIdentity:output:0*.
_input_shapes
:���������::20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:O K
'
_output_shapes
:���������
 
_user_specified_nameinputs
�	
�
B__inference_z_mean_layer_call_and_return_conditional_losses_141247

inputs"
matmul_readvariableop_resource#
biasadd_readvariableop_resource
identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOp�
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:*
dtype02
MatMul/ReadVariableOps
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2
MatMul�
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:*
dtype02
BiasAdd/ReadVariableOp�
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2	
BiasAdd�
IdentityIdentityBiasAdd:output:0^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*
T0*'
_output_shapes
:���������2

Identity"
identityIdentity:output:0*.
_input_shapes
:���������::20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:O K
'
_output_shapes
:���������
 
_user_specified_nameinputs
�

�
$__inference_signature_wrapper_142095
input_1
unknown
	unknown_0
	unknown_1
	unknown_2
	unknown_3
	unknown_4
	unknown_5
	unknown_6
	unknown_7
	unknown_8
	unknown_9

unknown_10

unknown_11

unknown_12

unknown_13

unknown_14

unknown_15

unknown_16
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinput_1unknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6	unknown_7	unknown_8	unknown_9
unknown_10
unknown_11
unknown_12
unknown_13
unknown_14
unknown_15
unknown_16*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*4
_read_only_resource_inputs
	
*-
config_proto

CPU

GPU 2J 8� **
f%R#
!__inference__wrapped_model_1411522
StatefulPartitionedCall�
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*'
_output_shapes
:���������2

Identity"
identityIdentity:output:0*n
_input_shapes]
[:���������::::::::::::::::::22
StatefulPartitionedCallStatefulPartitionedCall:P L
'
_output_shapes
:���������
!
_user_specified_name	input_1
�&
�
C__inference_encoder_layer_call_and_return_conditional_losses_141455

inputs
encoder_l2_141426
encoder_l2_141428
encoder_l3_141431
encoder_l3_141433
encoder_l4_141436
encoder_l4_141438
z_mean_141441
z_mean_141443
z_log_var_141446
z_log_var_141448
identity

identity_1

identity_2��"encoder_l2/StatefulPartitionedCall�"encoder_l3/StatefulPartitionedCall�"encoder_l4/StatefulPartitionedCall�"sampling_3/StatefulPartitionedCall�!z_log_var/StatefulPartitionedCall�z_mean/StatefulPartitionedCall�
"encoder_l2/StatefulPartitionedCallStatefulPartitionedCallinputsencoder_l2_141426encoder_l2_141428*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������
*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *O
fJRH
F__inference_encoder_l2_layer_call_and_return_conditional_losses_1411672$
"encoder_l2/StatefulPartitionedCall�
"encoder_l3/StatefulPartitionedCallStatefulPartitionedCall+encoder_l2/StatefulPartitionedCall:output:0encoder_l3_141431encoder_l3_141433*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *O
fJRH
F__inference_encoder_l3_layer_call_and_return_conditional_losses_1411942$
"encoder_l3/StatefulPartitionedCall�
"encoder_l4/StatefulPartitionedCallStatefulPartitionedCall+encoder_l3/StatefulPartitionedCall:output:0encoder_l4_141436encoder_l4_141438*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *O
fJRH
F__inference_encoder_l4_layer_call_and_return_conditional_losses_1412212$
"encoder_l4/StatefulPartitionedCall�
z_mean/StatefulPartitionedCallStatefulPartitionedCall+encoder_l4/StatefulPartitionedCall:output:0z_mean_141441z_mean_141443*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *K
fFRD
B__inference_z_mean_layer_call_and_return_conditional_losses_1412472 
z_mean/StatefulPartitionedCall�
!z_log_var/StatefulPartitionedCallStatefulPartitionedCall+encoder_l4/StatefulPartitionedCall:output:0z_log_var_141446z_log_var_141448*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *N
fIRG
E__inference_z_log_var_layer_call_and_return_conditional_losses_1412732#
!z_log_var/StatefulPartitionedCall�
"sampling_3/StatefulPartitionedCallStatefulPartitionedCall'z_mean/StatefulPartitionedCall:output:0*z_log_var/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8� *O
fJRH
F__inference_sampling_3_layer_call_and_return_conditional_losses_1413152$
"sampling_3/StatefulPartitionedCall�
IdentityIdentity'z_mean/StatefulPartitionedCall:output:0#^encoder_l2/StatefulPartitionedCall#^encoder_l3/StatefulPartitionedCall#^encoder_l4/StatefulPartitionedCall#^sampling_3/StatefulPartitionedCall"^z_log_var/StatefulPartitionedCall^z_mean/StatefulPartitionedCall*
T0*'
_output_shapes
:���������2

Identity�

Identity_1Identity*z_log_var/StatefulPartitionedCall:output:0#^encoder_l2/StatefulPartitionedCall#^encoder_l3/StatefulPartitionedCall#^encoder_l4/StatefulPartitionedCall#^sampling_3/StatefulPartitionedCall"^z_log_var/StatefulPartitionedCall^z_mean/StatefulPartitionedCall*
T0*'
_output_shapes
:���������2

Identity_1�

Identity_2Identity+sampling_3/StatefulPartitionedCall:output:0#^encoder_l2/StatefulPartitionedCall#^encoder_l3/StatefulPartitionedCall#^encoder_l4/StatefulPartitionedCall#^sampling_3/StatefulPartitionedCall"^z_log_var/StatefulPartitionedCall^z_mean/StatefulPartitionedCall*
T0*'
_output_shapes
:���������2

Identity_2"
identityIdentity:output:0"!

identity_1Identity_1:output:0"!

identity_2Identity_2:output:0*N
_input_shapes=
;:���������::::::::::2H
"encoder_l2/StatefulPartitionedCall"encoder_l2/StatefulPartitionedCall2H
"encoder_l3/StatefulPartitionedCall"encoder_l3/StatefulPartitionedCall2H
"encoder_l4/StatefulPartitionedCall"encoder_l4/StatefulPartitionedCall2H
"sampling_3/StatefulPartitionedCall"sampling_3/StatefulPartitionedCall2F
!z_log_var/StatefulPartitionedCall!z_log_var/StatefulPartitionedCall2@
z_mean/StatefulPartitionedCallz_mean/StatefulPartitionedCall:O K
'
_output_shapes
:���������
 
_user_specified_nameinputs
�R
�
C__inference_encoder_layer_call_and_return_conditional_losses_142747

inputs-
)encoder_l2_matmul_readvariableop_resource.
*encoder_l2_biasadd_readvariableop_resource-
)encoder_l3_matmul_readvariableop_resource.
*encoder_l3_biasadd_readvariableop_resource-
)encoder_l4_matmul_readvariableop_resource.
*encoder_l4_biasadd_readvariableop_resource)
%z_mean_matmul_readvariableop_resource*
&z_mean_biasadd_readvariableop_resource,
(z_log_var_matmul_readvariableop_resource-
)z_log_var_biasadd_readvariableop_resource
identity

identity_1

identity_2��!encoder_l2/BiasAdd/ReadVariableOp� encoder_l2/MatMul/ReadVariableOp�!encoder_l3/BiasAdd/ReadVariableOp� encoder_l3/MatMul/ReadVariableOp�!encoder_l4/BiasAdd/ReadVariableOp� encoder_l4/MatMul/ReadVariableOp� z_log_var/BiasAdd/ReadVariableOp�z_log_var/MatMul/ReadVariableOp�z_mean/BiasAdd/ReadVariableOp�z_mean/MatMul/ReadVariableOp�
 encoder_l2/MatMul/ReadVariableOpReadVariableOp)encoder_l2_matmul_readvariableop_resource*
_output_shapes

:
*
dtype02"
 encoder_l2/MatMul/ReadVariableOp�
encoder_l2/MatMulMatMulinputs(encoder_l2/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������
2
encoder_l2/MatMul�
!encoder_l2/BiasAdd/ReadVariableOpReadVariableOp*encoder_l2_biasadd_readvariableop_resource*
_output_shapes
:
*
dtype02#
!encoder_l2/BiasAdd/ReadVariableOp�
encoder_l2/BiasAddBiasAddencoder_l2/MatMul:product:0)encoder_l2/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������
2
encoder_l2/BiasAddy
encoder_l2/ReluReluencoder_l2/BiasAdd:output:0*
T0*'
_output_shapes
:���������
2
encoder_l2/Relu�
 encoder_l3/MatMul/ReadVariableOpReadVariableOp)encoder_l3_matmul_readvariableop_resource*
_output_shapes

:
*
dtype02"
 encoder_l3/MatMul/ReadVariableOp�
encoder_l3/MatMulMatMulencoder_l2/Relu:activations:0(encoder_l3/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2
encoder_l3/MatMul�
!encoder_l3/BiasAdd/ReadVariableOpReadVariableOp*encoder_l3_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02#
!encoder_l3/BiasAdd/ReadVariableOp�
encoder_l3/BiasAddBiasAddencoder_l3/MatMul:product:0)encoder_l3/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2
encoder_l3/BiasAddy
encoder_l3/ReluReluencoder_l3/BiasAdd:output:0*
T0*'
_output_shapes
:���������2
encoder_l3/Relu�
 encoder_l4/MatMul/ReadVariableOpReadVariableOp)encoder_l4_matmul_readvariableop_resource*
_output_shapes

:*
dtype02"
 encoder_l4/MatMul/ReadVariableOp�
encoder_l4/MatMulMatMulencoder_l3/Relu:activations:0(encoder_l4/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2
encoder_l4/MatMul�
!encoder_l4/BiasAdd/ReadVariableOpReadVariableOp*encoder_l4_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02#
!encoder_l4/BiasAdd/ReadVariableOp�
encoder_l4/BiasAddBiasAddencoder_l4/MatMul:product:0)encoder_l4/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2
encoder_l4/BiasAddy
encoder_l4/ReluReluencoder_l4/BiasAdd:output:0*
T0*'
_output_shapes
:���������2
encoder_l4/Relu�
z_mean/MatMul/ReadVariableOpReadVariableOp%z_mean_matmul_readvariableop_resource*
_output_shapes

:*
dtype02
z_mean/MatMul/ReadVariableOp�
z_mean/MatMulMatMulencoder_l4/Relu:activations:0$z_mean/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2
z_mean/MatMul�
z_mean/BiasAdd/ReadVariableOpReadVariableOp&z_mean_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02
z_mean/BiasAdd/ReadVariableOp�
z_mean/BiasAddBiasAddz_mean/MatMul:product:0%z_mean/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2
z_mean/BiasAdd�
z_log_var/MatMul/ReadVariableOpReadVariableOp(z_log_var_matmul_readvariableop_resource*
_output_shapes

:*
dtype02!
z_log_var/MatMul/ReadVariableOp�
z_log_var/MatMulMatMulencoder_l4/Relu:activations:0'z_log_var/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2
z_log_var/MatMul�
 z_log_var/BiasAdd/ReadVariableOpReadVariableOp)z_log_var_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02"
 z_log_var/BiasAdd/ReadVariableOp�
z_log_var/BiasAddBiasAddz_log_var/MatMul:product:0(z_log_var/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2
z_log_var/BiasAddk
sampling_3/ShapeShapez_mean/BiasAdd:output:0*
T0*
_output_shapes
:2
sampling_3/Shape�
sampling_3/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB: 2 
sampling_3/strided_slice/stack�
 sampling_3/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2"
 sampling_3/strided_slice/stack_1�
 sampling_3/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2"
 sampling_3/strided_slice/stack_2�
sampling_3/strided_sliceStridedSlicesampling_3/Shape:output:0'sampling_3/strided_slice/stack:output:0)sampling_3/strided_slice/stack_1:output:0)sampling_3/strided_slice/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2
sampling_3/strided_sliceo
sampling_3/Shape_1Shapez_mean/BiasAdd:output:0*
T0*
_output_shapes
:2
sampling_3/Shape_1�
 sampling_3/strided_slice_1/stackConst*
_output_shapes
:*
dtype0*
valueB:2"
 sampling_3/strided_slice_1/stack�
"sampling_3/strided_slice_1/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2$
"sampling_3/strided_slice_1/stack_1�
"sampling_3/strided_slice_1/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2$
"sampling_3/strided_slice_1/stack_2�
sampling_3/strided_slice_1StridedSlicesampling_3/Shape_1:output:0)sampling_3/strided_slice_1/stack:output:0+sampling_3/strided_slice_1/stack_1:output:0+sampling_3/strided_slice_1/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2
sampling_3/strided_slice_1�
sampling_3/random_normal/shapePack!sampling_3/strided_slice:output:0#sampling_3/strided_slice_1:output:0*
N*
T0*
_output_shapes
:2 
sampling_3/random_normal/shape�
sampling_3/random_normal/meanConst*
_output_shapes
: *
dtype0*
valueB
 *    2
sampling_3/random_normal/mean�
sampling_3/random_normal/stddevConst*
_output_shapes
: *
dtype0*
valueB
 *  �?2!
sampling_3/random_normal/stddev�
-sampling_3/random_normal/RandomStandardNormalRandomStandardNormal'sampling_3/random_normal/shape:output:0*
T0*0
_output_shapes
:������������������*
dtype0*
seed���)*
seed2̰�2/
-sampling_3/random_normal/RandomStandardNormal�
sampling_3/random_normal/mulMul6sampling_3/random_normal/RandomStandardNormal:output:0(sampling_3/random_normal/stddev:output:0*
T0*0
_output_shapes
:������������������2
sampling_3/random_normal/mul�
sampling_3/random_normalAdd sampling_3/random_normal/mul:z:0&sampling_3/random_normal/mean:output:0*
T0*0
_output_shapes
:������������������2
sampling_3/random_normali
sampling_3/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *   ?2
sampling_3/mul/x�
sampling_3/mulMulsampling_3/mul/x:output:0z_log_var/BiasAdd:output:0*
T0*'
_output_shapes
:���������2
sampling_3/mulm
sampling_3/ExpExpsampling_3/mul:z:0*
T0*'
_output_shapes
:���������2
sampling_3/Exp�
sampling_3/mul_1Mulsampling_3/Exp:y:0sampling_3/random_normal:z:0*
T0*'
_output_shapes
:���������2
sampling_3/mul_1�
sampling_3/addAddV2z_mean/BiasAdd:output:0sampling_3/mul_1:z:0*
T0*'
_output_shapes
:���������2
sampling_3/add�
IdentityIdentityz_mean/BiasAdd:output:0"^encoder_l2/BiasAdd/ReadVariableOp!^encoder_l2/MatMul/ReadVariableOp"^encoder_l3/BiasAdd/ReadVariableOp!^encoder_l3/MatMul/ReadVariableOp"^encoder_l4/BiasAdd/ReadVariableOp!^encoder_l4/MatMul/ReadVariableOp!^z_log_var/BiasAdd/ReadVariableOp ^z_log_var/MatMul/ReadVariableOp^z_mean/BiasAdd/ReadVariableOp^z_mean/MatMul/ReadVariableOp*
T0*'
_output_shapes
:���������2

Identity�

Identity_1Identityz_log_var/BiasAdd:output:0"^encoder_l2/BiasAdd/ReadVariableOp!^encoder_l2/MatMul/ReadVariableOp"^encoder_l3/BiasAdd/ReadVariableOp!^encoder_l3/MatMul/ReadVariableOp"^encoder_l4/BiasAdd/ReadVariableOp!^encoder_l4/MatMul/ReadVariableOp!^z_log_var/BiasAdd/ReadVariableOp ^z_log_var/MatMul/ReadVariableOp^z_mean/BiasAdd/ReadVariableOp^z_mean/MatMul/ReadVariableOp*
T0*'
_output_shapes
:���������2

Identity_1�

Identity_2Identitysampling_3/add:z:0"^encoder_l2/BiasAdd/ReadVariableOp!^encoder_l2/MatMul/ReadVariableOp"^encoder_l3/BiasAdd/ReadVariableOp!^encoder_l3/MatMul/ReadVariableOp"^encoder_l4/BiasAdd/ReadVariableOp!^encoder_l4/MatMul/ReadVariableOp!^z_log_var/BiasAdd/ReadVariableOp ^z_log_var/MatMul/ReadVariableOp^z_mean/BiasAdd/ReadVariableOp^z_mean/MatMul/ReadVariableOp*
T0*'
_output_shapes
:���������2

Identity_2"
identityIdentity:output:0"!

identity_1Identity_1:output:0"!

identity_2Identity_2:output:0*N
_input_shapes=
;:���������::::::::::2F
!encoder_l2/BiasAdd/ReadVariableOp!encoder_l2/BiasAdd/ReadVariableOp2D
 encoder_l2/MatMul/ReadVariableOp encoder_l2/MatMul/ReadVariableOp2F
!encoder_l3/BiasAdd/ReadVariableOp!encoder_l3/BiasAdd/ReadVariableOp2D
 encoder_l3/MatMul/ReadVariableOp encoder_l3/MatMul/ReadVariableOp2F
!encoder_l4/BiasAdd/ReadVariableOp!encoder_l4/BiasAdd/ReadVariableOp2D
 encoder_l4/MatMul/ReadVariableOp encoder_l4/MatMul/ReadVariableOp2D
 z_log_var/BiasAdd/ReadVariableOp z_log_var/BiasAdd/ReadVariableOp2B
z_log_var/MatMul/ReadVariableOpz_log_var/MatMul/ReadVariableOp2>
z_mean/BiasAdd/ReadVariableOpz_mean/BiasAdd/ReadVariableOp2<
z_mean/MatMul/ReadVariableOpz_mean/MatMul/ReadVariableOp:O K
'
_output_shapes
:���������
 
_user_specified_nameinputs
�s
�
__inference__traced_save_143327
file_prefix(
$savev2_adam_iter_read_readvariableop	*
&savev2_adam_beta_1_read_readvariableop*
&savev2_adam_beta_2_read_readvariableop)
%savev2_adam_decay_read_readvariableop1
-savev2_adam_learning_rate_read_readvariableop0
,savev2_encoder_l2_kernel_read_readvariableop.
*savev2_encoder_l2_bias_read_readvariableop0
,savev2_encoder_l3_kernel_read_readvariableop.
*savev2_encoder_l3_bias_read_readvariableop0
,savev2_encoder_l4_kernel_read_readvariableop.
*savev2_encoder_l4_bias_read_readvariableop,
(savev2_z_mean_kernel_read_readvariableop*
&savev2_z_mean_bias_read_readvariableop/
+savev2_z_log_var_kernel_read_readvariableop-
)savev2_z_log_var_bias_read_readvariableop0
,savev2_decoder_l2_kernel_read_readvariableop.
*savev2_decoder_l2_bias_read_readvariableop0
,savev2_decoder_l3_kernel_read_readvariableop.
*savev2_decoder_l3_bias_read_readvariableop0
,savev2_decoder_l4_kernel_read_readvariableop.
*savev2_decoder_l4_bias_read_readvariableop2
.savev2_output_layer_kernel_read_readvariableop0
,savev2_output_layer_bias_read_readvariableop$
 savev2_total_read_readvariableop$
 savev2_count_read_readvariableop7
3savev2_adam_encoder_l2_kernel_m_read_readvariableop5
1savev2_adam_encoder_l2_bias_m_read_readvariableop7
3savev2_adam_encoder_l3_kernel_m_read_readvariableop5
1savev2_adam_encoder_l3_bias_m_read_readvariableop7
3savev2_adam_encoder_l4_kernel_m_read_readvariableop5
1savev2_adam_encoder_l4_bias_m_read_readvariableop3
/savev2_adam_z_mean_kernel_m_read_readvariableop1
-savev2_adam_z_mean_bias_m_read_readvariableop6
2savev2_adam_z_log_var_kernel_m_read_readvariableop4
0savev2_adam_z_log_var_bias_m_read_readvariableop7
3savev2_adam_decoder_l2_kernel_m_read_readvariableop5
1savev2_adam_decoder_l2_bias_m_read_readvariableop7
3savev2_adam_decoder_l3_kernel_m_read_readvariableop5
1savev2_adam_decoder_l3_bias_m_read_readvariableop7
3savev2_adam_decoder_l4_kernel_m_read_readvariableop5
1savev2_adam_decoder_l4_bias_m_read_readvariableop9
5savev2_adam_output_layer_kernel_m_read_readvariableop7
3savev2_adam_output_layer_bias_m_read_readvariableop7
3savev2_adam_encoder_l2_kernel_v_read_readvariableop5
1savev2_adam_encoder_l2_bias_v_read_readvariableop7
3savev2_adam_encoder_l3_kernel_v_read_readvariableop5
1savev2_adam_encoder_l3_bias_v_read_readvariableop7
3savev2_adam_encoder_l4_kernel_v_read_readvariableop5
1savev2_adam_encoder_l4_bias_v_read_readvariableop3
/savev2_adam_z_mean_kernel_v_read_readvariableop1
-savev2_adam_z_mean_bias_v_read_readvariableop6
2savev2_adam_z_log_var_kernel_v_read_readvariableop4
0savev2_adam_z_log_var_bias_v_read_readvariableop7
3savev2_adam_decoder_l2_kernel_v_read_readvariableop5
1savev2_adam_decoder_l2_bias_v_read_readvariableop7
3savev2_adam_decoder_l3_kernel_v_read_readvariableop5
1savev2_adam_decoder_l3_bias_v_read_readvariableop7
3savev2_adam_decoder_l4_kernel_v_read_readvariableop5
1savev2_adam_decoder_l4_bias_v_read_readvariableop9
5savev2_adam_output_layer_kernel_v_read_readvariableop7
3savev2_adam_output_layer_bias_v_read_readvariableop
savev2_const

identity_1��MergeV2Checkpoints�
StaticRegexFullMatchStaticRegexFullMatchfile_prefix"/device:CPU:**
_output_shapes
: *
pattern
^s3://.*2
StaticRegexFullMatchc
ConstConst"/device:CPU:**
_output_shapes
: *
dtype0*
valueB B.part2
Constl
Const_1Const"/device:CPU:**
_output_shapes
: *
dtype0*
valueB B
_temp/part2	
Const_1�
SelectSelectStaticRegexFullMatch:output:0Const:output:0Const_1:output:0"/device:CPU:**
T0*
_output_shapes
: 2
Selectt

StringJoin
StringJoinfile_prefixSelect:output:0"/device:CPU:**
N*
_output_shapes
: 2

StringJoinZ

num_shardsConst*
_output_shapes
: *
dtype0*
value	B :2

num_shards
ShardedFilename/shardConst"/device:CPU:0*
_output_shapes
: *
dtype0*
value	B : 2
ShardedFilename/shard�
ShardedFilenameShardedFilenameStringJoin:output:0ShardedFilename/shard:output:0num_shards:output:0"/device:CPU:0*
_output_shapes
: 2
ShardedFilename�
SaveV2/tensor_namesConst"/device:CPU:0*
_output_shapes
:>*
dtype0*�
value�B�>B)optimizer/iter/.ATTRIBUTES/VARIABLE_VALUEB+optimizer/beta_1/.ATTRIBUTES/VARIABLE_VALUEB+optimizer/beta_2/.ATTRIBUTES/VARIABLE_VALUEB*optimizer/decay/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/learning_rate/.ATTRIBUTES/VARIABLE_VALUEB&variables/0/.ATTRIBUTES/VARIABLE_VALUEB&variables/1/.ATTRIBUTES/VARIABLE_VALUEB&variables/2/.ATTRIBUTES/VARIABLE_VALUEB&variables/3/.ATTRIBUTES/VARIABLE_VALUEB&variables/4/.ATTRIBUTES/VARIABLE_VALUEB&variables/5/.ATTRIBUTES/VARIABLE_VALUEB&variables/6/.ATTRIBUTES/VARIABLE_VALUEB&variables/7/.ATTRIBUTES/VARIABLE_VALUEB&variables/8/.ATTRIBUTES/VARIABLE_VALUEB&variables/9/.ATTRIBUTES/VARIABLE_VALUEB'variables/10/.ATTRIBUTES/VARIABLE_VALUEB'variables/11/.ATTRIBUTES/VARIABLE_VALUEB'variables/12/.ATTRIBUTES/VARIABLE_VALUEB'variables/13/.ATTRIBUTES/VARIABLE_VALUEB'variables/14/.ATTRIBUTES/VARIABLE_VALUEB'variables/15/.ATTRIBUTES/VARIABLE_VALUEB'variables/16/.ATTRIBUTES/VARIABLE_VALUEB'variables/17/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/0/total/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/0/count/.ATTRIBUTES/VARIABLE_VALUEBBvariables/0/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBBvariables/1/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBBvariables/2/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBBvariables/3/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBBvariables/4/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBBvariables/5/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBBvariables/6/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBBvariables/7/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBBvariables/8/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBBvariables/9/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBCvariables/10/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBCvariables/11/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBCvariables/12/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBCvariables/13/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBCvariables/14/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBCvariables/15/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBCvariables/16/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBCvariables/17/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBBvariables/0/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBBvariables/1/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBBvariables/2/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBBvariables/3/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBBvariables/4/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBBvariables/5/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBBvariables/6/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBBvariables/7/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBBvariables/8/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBBvariables/9/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBCvariables/10/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBCvariables/11/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBCvariables/12/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBCvariables/13/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBCvariables/14/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBCvariables/15/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBCvariables/16/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBCvariables/17/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEB_CHECKPOINTABLE_OBJECT_GRAPH2
SaveV2/tensor_names�
SaveV2/shape_and_slicesConst"/device:CPU:0*
_output_shapes
:>*
dtype0*�
value�B�>B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B 2
SaveV2/shape_and_slices�
SaveV2SaveV2ShardedFilename:filename:0SaveV2/tensor_names:output:0 SaveV2/shape_and_slices:output:0$savev2_adam_iter_read_readvariableop&savev2_adam_beta_1_read_readvariableop&savev2_adam_beta_2_read_readvariableop%savev2_adam_decay_read_readvariableop-savev2_adam_learning_rate_read_readvariableop,savev2_encoder_l2_kernel_read_readvariableop*savev2_encoder_l2_bias_read_readvariableop,savev2_encoder_l3_kernel_read_readvariableop*savev2_encoder_l3_bias_read_readvariableop,savev2_encoder_l4_kernel_read_readvariableop*savev2_encoder_l4_bias_read_readvariableop(savev2_z_mean_kernel_read_readvariableop&savev2_z_mean_bias_read_readvariableop+savev2_z_log_var_kernel_read_readvariableop)savev2_z_log_var_bias_read_readvariableop,savev2_decoder_l2_kernel_read_readvariableop*savev2_decoder_l2_bias_read_readvariableop,savev2_decoder_l3_kernel_read_readvariableop*savev2_decoder_l3_bias_read_readvariableop,savev2_decoder_l4_kernel_read_readvariableop*savev2_decoder_l4_bias_read_readvariableop.savev2_output_layer_kernel_read_readvariableop,savev2_output_layer_bias_read_readvariableop savev2_total_read_readvariableop savev2_count_read_readvariableop3savev2_adam_encoder_l2_kernel_m_read_readvariableop1savev2_adam_encoder_l2_bias_m_read_readvariableop3savev2_adam_encoder_l3_kernel_m_read_readvariableop1savev2_adam_encoder_l3_bias_m_read_readvariableop3savev2_adam_encoder_l4_kernel_m_read_readvariableop1savev2_adam_encoder_l4_bias_m_read_readvariableop/savev2_adam_z_mean_kernel_m_read_readvariableop-savev2_adam_z_mean_bias_m_read_readvariableop2savev2_adam_z_log_var_kernel_m_read_readvariableop0savev2_adam_z_log_var_bias_m_read_readvariableop3savev2_adam_decoder_l2_kernel_m_read_readvariableop1savev2_adam_decoder_l2_bias_m_read_readvariableop3savev2_adam_decoder_l3_kernel_m_read_readvariableop1savev2_adam_decoder_l3_bias_m_read_readvariableop3savev2_adam_decoder_l4_kernel_m_read_readvariableop1savev2_adam_decoder_l4_bias_m_read_readvariableop5savev2_adam_output_layer_kernel_m_read_readvariableop3savev2_adam_output_layer_bias_m_read_readvariableop3savev2_adam_encoder_l2_kernel_v_read_readvariableop1savev2_adam_encoder_l2_bias_v_read_readvariableop3savev2_adam_encoder_l3_kernel_v_read_readvariableop1savev2_adam_encoder_l3_bias_v_read_readvariableop3savev2_adam_encoder_l4_kernel_v_read_readvariableop1savev2_adam_encoder_l4_bias_v_read_readvariableop/savev2_adam_z_mean_kernel_v_read_readvariableop-savev2_adam_z_mean_bias_v_read_readvariableop2savev2_adam_z_log_var_kernel_v_read_readvariableop0savev2_adam_z_log_var_bias_v_read_readvariableop3savev2_adam_decoder_l2_kernel_v_read_readvariableop1savev2_adam_decoder_l2_bias_v_read_readvariableop3savev2_adam_decoder_l3_kernel_v_read_readvariableop1savev2_adam_decoder_l3_bias_v_read_readvariableop3savev2_adam_decoder_l4_kernel_v_read_readvariableop1savev2_adam_decoder_l4_bias_v_read_readvariableop5savev2_adam_output_layer_kernel_v_read_readvariableop3savev2_adam_output_layer_bias_v_read_readvariableopsavev2_const"/device:CPU:0*
_output_shapes
 *L
dtypesB
@2>	2
SaveV2�
&MergeV2Checkpoints/checkpoint_prefixesPackShardedFilename:filename:0^SaveV2"/device:CPU:0*
N*
T0*
_output_shapes
:2(
&MergeV2Checkpoints/checkpoint_prefixes�
MergeV2CheckpointsMergeV2Checkpoints/MergeV2Checkpoints/checkpoint_prefixes:output:0file_prefix"/device:CPU:0*
_output_shapes
 2
MergeV2Checkpointsr
IdentityIdentityfile_prefix^MergeV2Checkpoints"/device:CPU:0*
T0*
_output_shapes
: 2

Identitym

Identity_1IdentityIdentity:output:0^MergeV2Checkpoints*
T0*
_output_shapes
: 2

Identity_1"!

identity_1Identity_1:output:0*�
_input_shapes�
�: : : : : : :
:
:
::::::::::::
:
:
:: : :
:
:
::::::::::::
:
:
::
:
:
::::::::::::
:
:
:: 2(
MergeV2CheckpointsMergeV2Checkpoints:C ?

_output_shapes
: 
%
_user_specified_namefile_prefix:

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :$ 

_output_shapes

:
: 

_output_shapes
:
:$ 

_output_shapes

:
: 	

_output_shapes
::$
 

_output_shapes

:: 

_output_shapes
::$ 

_output_shapes

:: 

_output_shapes
::$ 

_output_shapes

:: 

_output_shapes
::$ 

_output_shapes

:: 

_output_shapes
::$ 

_output_shapes

:: 

_output_shapes
::$ 

_output_shapes

:
: 

_output_shapes
:
:$ 

_output_shapes

:
: 

_output_shapes
::

_output_shapes
: :

_output_shapes
: :$ 

_output_shapes

:
: 

_output_shapes
:
:$ 

_output_shapes

:
: 

_output_shapes
::$ 

_output_shapes

:: 

_output_shapes
::$  

_output_shapes

:: !

_output_shapes
::$" 

_output_shapes

:: #

_output_shapes
::$$ 

_output_shapes

:: %

_output_shapes
::$& 

_output_shapes

:: '

_output_shapes
::$( 

_output_shapes

:
: )

_output_shapes
:
:$* 

_output_shapes

:
: +

_output_shapes
::$, 

_output_shapes

:
: -

_output_shapes
:
:$. 

_output_shapes

:
: /

_output_shapes
::$0 

_output_shapes

:: 1

_output_shapes
::$2 

_output_shapes

:: 3

_output_shapes
::$4 

_output_shapes

:: 5

_output_shapes
::$6 

_output_shapes

:: 7

_output_shapes
::$8 

_output_shapes

:: 9

_output_shapes
::$: 

_output_shapes

:
: ;

_output_shapes
:
:$< 

_output_shapes

:
: =

_output_shapes
::>

_output_shapes
: 
��
�
B__inference_cvae_3_layer_call_and_return_conditional_losses_142545

inputs5
1encoder_encoder_l2_matmul_readvariableop_resource6
2encoder_encoder_l2_biasadd_readvariableop_resource5
1encoder_encoder_l3_matmul_readvariableop_resource6
2encoder_encoder_l3_biasadd_readvariableop_resource5
1encoder_encoder_l4_matmul_readvariableop_resource6
2encoder_encoder_l4_biasadd_readvariableop_resource1
-encoder_z_mean_matmul_readvariableop_resource2
.encoder_z_mean_biasadd_readvariableop_resource4
0encoder_z_log_var_matmul_readvariableop_resource5
1encoder_z_log_var_biasadd_readvariableop_resource5
1decoder_decoder_l2_matmul_readvariableop_resource6
2decoder_decoder_l2_biasadd_readvariableop_resource5
1decoder_decoder_l3_matmul_readvariableop_resource6
2decoder_decoder_l3_biasadd_readvariableop_resource5
1decoder_decoder_l4_matmul_readvariableop_resource6
2decoder_decoder_l4_biasadd_readvariableop_resource7
3decoder_output_layer_matmul_readvariableop_resource8
4decoder_output_layer_biasadd_readvariableop_resource
identity��)decoder/decoder_l2/BiasAdd/ReadVariableOp�(decoder/decoder_l2/MatMul/ReadVariableOp�)decoder/decoder_l3/BiasAdd/ReadVariableOp�(decoder/decoder_l3/MatMul/ReadVariableOp�)decoder/decoder_l4/BiasAdd/ReadVariableOp�(decoder/decoder_l4/MatMul/ReadVariableOp�+decoder/output_layer/BiasAdd/ReadVariableOp�*decoder/output_layer/MatMul/ReadVariableOp�)encoder/encoder_l2/BiasAdd/ReadVariableOp�(encoder/encoder_l2/MatMul/ReadVariableOp�)encoder/encoder_l3/BiasAdd/ReadVariableOp�(encoder/encoder_l3/MatMul/ReadVariableOp�)encoder/encoder_l4/BiasAdd/ReadVariableOp�(encoder/encoder_l4/MatMul/ReadVariableOp�(encoder/z_log_var/BiasAdd/ReadVariableOp�'encoder/z_log_var/MatMul/ReadVariableOp�%encoder/z_mean/BiasAdd/ReadVariableOp�$encoder/z_mean/MatMul/ReadVariableOp�
(encoder/encoder_l2/MatMul/ReadVariableOpReadVariableOp1encoder_encoder_l2_matmul_readvariableop_resource*
_output_shapes

:
*
dtype02*
(encoder/encoder_l2/MatMul/ReadVariableOp�
encoder/encoder_l2/MatMulMatMulinputs0encoder/encoder_l2/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������
2
encoder/encoder_l2/MatMul�
)encoder/encoder_l2/BiasAdd/ReadVariableOpReadVariableOp2encoder_encoder_l2_biasadd_readvariableop_resource*
_output_shapes
:
*
dtype02+
)encoder/encoder_l2/BiasAdd/ReadVariableOp�
encoder/encoder_l2/BiasAddBiasAdd#encoder/encoder_l2/MatMul:product:01encoder/encoder_l2/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������
2
encoder/encoder_l2/BiasAdd�
encoder/encoder_l2/ReluRelu#encoder/encoder_l2/BiasAdd:output:0*
T0*'
_output_shapes
:���������
2
encoder/encoder_l2/Relu�
(encoder/encoder_l3/MatMul/ReadVariableOpReadVariableOp1encoder_encoder_l3_matmul_readvariableop_resource*
_output_shapes

:
*
dtype02*
(encoder/encoder_l3/MatMul/ReadVariableOp�
encoder/encoder_l3/MatMulMatMul%encoder/encoder_l2/Relu:activations:00encoder/encoder_l3/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2
encoder/encoder_l3/MatMul�
)encoder/encoder_l3/BiasAdd/ReadVariableOpReadVariableOp2encoder_encoder_l3_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02+
)encoder/encoder_l3/BiasAdd/ReadVariableOp�
encoder/encoder_l3/BiasAddBiasAdd#encoder/encoder_l3/MatMul:product:01encoder/encoder_l3/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2
encoder/encoder_l3/BiasAdd�
encoder/encoder_l3/ReluRelu#encoder/encoder_l3/BiasAdd:output:0*
T0*'
_output_shapes
:���������2
encoder/encoder_l3/Relu�
(encoder/encoder_l4/MatMul/ReadVariableOpReadVariableOp1encoder_encoder_l4_matmul_readvariableop_resource*
_output_shapes

:*
dtype02*
(encoder/encoder_l4/MatMul/ReadVariableOp�
encoder/encoder_l4/MatMulMatMul%encoder/encoder_l3/Relu:activations:00encoder/encoder_l4/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2
encoder/encoder_l4/MatMul�
)encoder/encoder_l4/BiasAdd/ReadVariableOpReadVariableOp2encoder_encoder_l4_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02+
)encoder/encoder_l4/BiasAdd/ReadVariableOp�
encoder/encoder_l4/BiasAddBiasAdd#encoder/encoder_l4/MatMul:product:01encoder/encoder_l4/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2
encoder/encoder_l4/BiasAdd�
encoder/encoder_l4/ReluRelu#encoder/encoder_l4/BiasAdd:output:0*
T0*'
_output_shapes
:���������2
encoder/encoder_l4/Relu�
$encoder/z_mean/MatMul/ReadVariableOpReadVariableOp-encoder_z_mean_matmul_readvariableop_resource*
_output_shapes

:*
dtype02&
$encoder/z_mean/MatMul/ReadVariableOp�
encoder/z_mean/MatMulMatMul%encoder/encoder_l4/Relu:activations:0,encoder/z_mean/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2
encoder/z_mean/MatMul�
%encoder/z_mean/BiasAdd/ReadVariableOpReadVariableOp.encoder_z_mean_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02'
%encoder/z_mean/BiasAdd/ReadVariableOp�
encoder/z_mean/BiasAddBiasAddencoder/z_mean/MatMul:product:0-encoder/z_mean/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2
encoder/z_mean/BiasAdd�
'encoder/z_log_var/MatMul/ReadVariableOpReadVariableOp0encoder_z_log_var_matmul_readvariableop_resource*
_output_shapes

:*
dtype02)
'encoder/z_log_var/MatMul/ReadVariableOp�
encoder/z_log_var/MatMulMatMul%encoder/encoder_l4/Relu:activations:0/encoder/z_log_var/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2
encoder/z_log_var/MatMul�
(encoder/z_log_var/BiasAdd/ReadVariableOpReadVariableOp1encoder_z_log_var_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02*
(encoder/z_log_var/BiasAdd/ReadVariableOp�
encoder/z_log_var/BiasAddBiasAdd"encoder/z_log_var/MatMul:product:00encoder/z_log_var/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2
encoder/z_log_var/BiasAdd�
encoder/sampling_3/ShapeShapeencoder/z_mean/BiasAdd:output:0*
T0*
_output_shapes
:2
encoder/sampling_3/Shape�
&encoder/sampling_3/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB: 2(
&encoder/sampling_3/strided_slice/stack�
(encoder/sampling_3/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2*
(encoder/sampling_3/strided_slice/stack_1�
(encoder/sampling_3/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2*
(encoder/sampling_3/strided_slice/stack_2�
 encoder/sampling_3/strided_sliceStridedSlice!encoder/sampling_3/Shape:output:0/encoder/sampling_3/strided_slice/stack:output:01encoder/sampling_3/strided_slice/stack_1:output:01encoder/sampling_3/strided_slice/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2"
 encoder/sampling_3/strided_slice�
encoder/sampling_3/Shape_1Shapeencoder/z_mean/BiasAdd:output:0*
T0*
_output_shapes
:2
encoder/sampling_3/Shape_1�
(encoder/sampling_3/strided_slice_1/stackConst*
_output_shapes
:*
dtype0*
valueB:2*
(encoder/sampling_3/strided_slice_1/stack�
*encoder/sampling_3/strided_slice_1/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2,
*encoder/sampling_3/strided_slice_1/stack_1�
*encoder/sampling_3/strided_slice_1/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2,
*encoder/sampling_3/strided_slice_1/stack_2�
"encoder/sampling_3/strided_slice_1StridedSlice#encoder/sampling_3/Shape_1:output:01encoder/sampling_3/strided_slice_1/stack:output:03encoder/sampling_3/strided_slice_1/stack_1:output:03encoder/sampling_3/strided_slice_1/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2$
"encoder/sampling_3/strided_slice_1�
&encoder/sampling_3/random_normal/shapePack)encoder/sampling_3/strided_slice:output:0+encoder/sampling_3/strided_slice_1:output:0*
N*
T0*
_output_shapes
:2(
&encoder/sampling_3/random_normal/shape�
%encoder/sampling_3/random_normal/meanConst*
_output_shapes
: *
dtype0*
valueB
 *    2'
%encoder/sampling_3/random_normal/mean�
'encoder/sampling_3/random_normal/stddevConst*
_output_shapes
: *
dtype0*
valueB
 *  �?2)
'encoder/sampling_3/random_normal/stddev�
5encoder/sampling_3/random_normal/RandomStandardNormalRandomStandardNormal/encoder/sampling_3/random_normal/shape:output:0*
T0*0
_output_shapes
:������������������*
dtype0*
seed���)*
seed2�_27
5encoder/sampling_3/random_normal/RandomStandardNormal�
$encoder/sampling_3/random_normal/mulMul>encoder/sampling_3/random_normal/RandomStandardNormal:output:00encoder/sampling_3/random_normal/stddev:output:0*
T0*0
_output_shapes
:������������������2&
$encoder/sampling_3/random_normal/mul�
 encoder/sampling_3/random_normalAdd(encoder/sampling_3/random_normal/mul:z:0.encoder/sampling_3/random_normal/mean:output:0*
T0*0
_output_shapes
:������������������2"
 encoder/sampling_3/random_normaly
encoder/sampling_3/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *   ?2
encoder/sampling_3/mul/x�
encoder/sampling_3/mulMul!encoder/sampling_3/mul/x:output:0"encoder/z_log_var/BiasAdd:output:0*
T0*'
_output_shapes
:���������2
encoder/sampling_3/mul�
encoder/sampling_3/ExpExpencoder/sampling_3/mul:z:0*
T0*'
_output_shapes
:���������2
encoder/sampling_3/Exp�
encoder/sampling_3/mul_1Mulencoder/sampling_3/Exp:y:0$encoder/sampling_3/random_normal:z:0*
T0*'
_output_shapes
:���������2
encoder/sampling_3/mul_1�
encoder/sampling_3/addAddV2encoder/z_mean/BiasAdd:output:0encoder/sampling_3/mul_1:z:0*
T0*'
_output_shapes
:���������2
encoder/sampling_3/add{
strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB"    ����2
strided_slice/stack
strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB"        2
strided_slice/stack_1
strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      2
strided_slice/stack_2�
strided_sliceStridedSliceinputsstrided_slice/stack:output:0strided_slice/stack_1:output:0strided_slice/stack_2:output:0*
Index0*
T0*'
_output_shapes
:���������*

begin_mask*
end_mask2
strided_slicet
concatenate/concat/axisConst*
_output_shapes
: *
dtype0*
value	B :2
concatenate/concat/axis�
concatenate/concatConcatV2encoder/sampling_3/add:z:0strided_slice:output:0 concatenate/concat/axis:output:0*
N*
T0*'
_output_shapes
:���������2
concatenate/concat�
(decoder/decoder_l2/MatMul/ReadVariableOpReadVariableOp1decoder_decoder_l2_matmul_readvariableop_resource*
_output_shapes

:*
dtype02*
(decoder/decoder_l2/MatMul/ReadVariableOp�
decoder/decoder_l2/MatMulMatMulconcatenate/concat:output:00decoder/decoder_l2/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2
decoder/decoder_l2/MatMul�
)decoder/decoder_l2/BiasAdd/ReadVariableOpReadVariableOp2decoder_decoder_l2_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02+
)decoder/decoder_l2/BiasAdd/ReadVariableOp�
decoder/decoder_l2/BiasAddBiasAdd#decoder/decoder_l2/MatMul:product:01decoder/decoder_l2/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2
decoder/decoder_l2/BiasAdd�
decoder/decoder_l2/ReluRelu#decoder/decoder_l2/BiasAdd:output:0*
T0*'
_output_shapes
:���������2
decoder/decoder_l2/Relu�
(decoder/decoder_l3/MatMul/ReadVariableOpReadVariableOp1decoder_decoder_l3_matmul_readvariableop_resource*
_output_shapes

:*
dtype02*
(decoder/decoder_l3/MatMul/ReadVariableOp�
decoder/decoder_l3/MatMulMatMul%decoder/decoder_l2/Relu:activations:00decoder/decoder_l3/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2
decoder/decoder_l3/MatMul�
)decoder/decoder_l3/BiasAdd/ReadVariableOpReadVariableOp2decoder_decoder_l3_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02+
)decoder/decoder_l3/BiasAdd/ReadVariableOp�
decoder/decoder_l3/BiasAddBiasAdd#decoder/decoder_l3/MatMul:product:01decoder/decoder_l3/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2
decoder/decoder_l3/BiasAdd�
decoder/decoder_l3/ReluRelu#decoder/decoder_l3/BiasAdd:output:0*
T0*'
_output_shapes
:���������2
decoder/decoder_l3/Relu�
(decoder/decoder_l4/MatMul/ReadVariableOpReadVariableOp1decoder_decoder_l4_matmul_readvariableop_resource*
_output_shapes

:
*
dtype02*
(decoder/decoder_l4/MatMul/ReadVariableOp�
decoder/decoder_l4/MatMulMatMul%decoder/decoder_l3/Relu:activations:00decoder/decoder_l4/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������
2
decoder/decoder_l4/MatMul�
)decoder/decoder_l4/BiasAdd/ReadVariableOpReadVariableOp2decoder_decoder_l4_biasadd_readvariableop_resource*
_output_shapes
:
*
dtype02+
)decoder/decoder_l4/BiasAdd/ReadVariableOp�
decoder/decoder_l4/BiasAddBiasAdd#decoder/decoder_l4/MatMul:product:01decoder/decoder_l4/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������
2
decoder/decoder_l4/BiasAdd�
decoder/decoder_l4/ReluRelu#decoder/decoder_l4/BiasAdd:output:0*
T0*'
_output_shapes
:���������
2
decoder/decoder_l4/Relu�
*decoder/output_layer/MatMul/ReadVariableOpReadVariableOp3decoder_output_layer_matmul_readvariableop_resource*
_output_shapes

:
*
dtype02,
*decoder/output_layer/MatMul/ReadVariableOp�
decoder/output_layer/MatMulMatMul%decoder/decoder_l4/Relu:activations:02decoder/output_layer/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2
decoder/output_layer/MatMul�
+decoder/output_layer/BiasAdd/ReadVariableOpReadVariableOp4decoder_output_layer_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02-
+decoder/output_layer/BiasAdd/ReadVariableOp�
decoder/output_layer/BiasAddBiasAdd%decoder/output_layer/MatMul:product:03decoder/output_layer/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2
decoder/output_layer/BiasAdd�
decoder/output_layer/ReluRelu%decoder/output_layer/BiasAdd:output:0*
T0*'
_output_shapes
:���������2
decoder/output_layer/Relu�
IdentityIdentity'decoder/output_layer/Relu:activations:0*^decoder/decoder_l2/BiasAdd/ReadVariableOp)^decoder/decoder_l2/MatMul/ReadVariableOp*^decoder/decoder_l3/BiasAdd/ReadVariableOp)^decoder/decoder_l3/MatMul/ReadVariableOp*^decoder/decoder_l4/BiasAdd/ReadVariableOp)^decoder/decoder_l4/MatMul/ReadVariableOp,^decoder/output_layer/BiasAdd/ReadVariableOp+^decoder/output_layer/MatMul/ReadVariableOp*^encoder/encoder_l2/BiasAdd/ReadVariableOp)^encoder/encoder_l2/MatMul/ReadVariableOp*^encoder/encoder_l3/BiasAdd/ReadVariableOp)^encoder/encoder_l3/MatMul/ReadVariableOp*^encoder/encoder_l4/BiasAdd/ReadVariableOp)^encoder/encoder_l4/MatMul/ReadVariableOp)^encoder/z_log_var/BiasAdd/ReadVariableOp(^encoder/z_log_var/MatMul/ReadVariableOp&^encoder/z_mean/BiasAdd/ReadVariableOp%^encoder/z_mean/MatMul/ReadVariableOp*
T0*'
_output_shapes
:���������2

Identity"
identityIdentity:output:0*n
_input_shapes]
[:���������::::::::::::::::::2V
)decoder/decoder_l2/BiasAdd/ReadVariableOp)decoder/decoder_l2/BiasAdd/ReadVariableOp2T
(decoder/decoder_l2/MatMul/ReadVariableOp(decoder/decoder_l2/MatMul/ReadVariableOp2V
)decoder/decoder_l3/BiasAdd/ReadVariableOp)decoder/decoder_l3/BiasAdd/ReadVariableOp2T
(decoder/decoder_l3/MatMul/ReadVariableOp(decoder/decoder_l3/MatMul/ReadVariableOp2V
)decoder/decoder_l4/BiasAdd/ReadVariableOp)decoder/decoder_l4/BiasAdd/ReadVariableOp2T
(decoder/decoder_l4/MatMul/ReadVariableOp(decoder/decoder_l4/MatMul/ReadVariableOp2Z
+decoder/output_layer/BiasAdd/ReadVariableOp+decoder/output_layer/BiasAdd/ReadVariableOp2X
*decoder/output_layer/MatMul/ReadVariableOp*decoder/output_layer/MatMul/ReadVariableOp2V
)encoder/encoder_l2/BiasAdd/ReadVariableOp)encoder/encoder_l2/BiasAdd/ReadVariableOp2T
(encoder/encoder_l2/MatMul/ReadVariableOp(encoder/encoder_l2/MatMul/ReadVariableOp2V
)encoder/encoder_l3/BiasAdd/ReadVariableOp)encoder/encoder_l3/BiasAdd/ReadVariableOp2T
(encoder/encoder_l3/MatMul/ReadVariableOp(encoder/encoder_l3/MatMul/ReadVariableOp2V
)encoder/encoder_l4/BiasAdd/ReadVariableOp)encoder/encoder_l4/BiasAdd/ReadVariableOp2T
(encoder/encoder_l4/MatMul/ReadVariableOp(encoder/encoder_l4/MatMul/ReadVariableOp2T
(encoder/z_log_var/BiasAdd/ReadVariableOp(encoder/z_log_var/BiasAdd/ReadVariableOp2R
'encoder/z_log_var/MatMul/ReadVariableOp'encoder/z_log_var/MatMul/ReadVariableOp2N
%encoder/z_mean/BiasAdd/ReadVariableOp%encoder/z_mean/BiasAdd/ReadVariableOp2L
$encoder/z_mean/MatMul/ReadVariableOp$encoder/z_mean/MatMul/ReadVariableOp:O K
'
_output_shapes
:���������
 
_user_specified_nameinputs
�

*__inference_z_log_var_layer_call_fn_143009

inputs
unknown
	unknown_0
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *N
fIRG
E__inference_z_log_var_layer_call_and_return_conditional_losses_1412732
StatefulPartitionedCall�
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*'
_output_shapes
:���������2

Identity"
identityIdentity:output:0*.
_input_shapes
:���������::22
StatefulPartitionedCallStatefulPartitionedCall:O K
'
_output_shapes
:���������
 
_user_specified_nameinputs
�

�
'__inference_cvae_3_layer_call_fn_142627

inputs
unknown
	unknown_0
	unknown_1
	unknown_2
	unknown_3
	unknown_4
	unknown_5
	unknown_6
	unknown_7
	unknown_8
	unknown_9

unknown_10

unknown_11

unknown_12

unknown_13

unknown_14

unknown_15

unknown_16
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6	unknown_7	unknown_8	unknown_9
unknown_10
unknown_11
unknown_12
unknown_13
unknown_14
unknown_15
unknown_16*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*4
_read_only_resource_inputs
	
*-
config_proto

CPU

GPU 2J 8� *K
fFRD
B__inference_cvae_3_layer_call_and_return_conditional_losses_1419642
StatefulPartitionedCall�
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*'
_output_shapes
:���������2

Identity"
identityIdentity:output:0*n
_input_shapes]
[:���������::::::::::::::::::22
StatefulPartitionedCallStatefulPartitionedCall:O K
'
_output_shapes
:���������
 
_user_specified_nameinputs
�

�
'__inference_cvae_3_layer_call_fn_142586

inputs
unknown
	unknown_0
	unknown_1
	unknown_2
	unknown_3
	unknown_4
	unknown_5
	unknown_6
	unknown_7
	unknown_8
	unknown_9

unknown_10

unknown_11

unknown_12

unknown_13

unknown_14

unknown_15

unknown_16
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6	unknown_7	unknown_8	unknown_9
unknown_10
unknown_11
unknown_12
unknown_13
unknown_14
unknown_15
unknown_16*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*4
_read_only_resource_inputs
	
*-
config_proto

CPU

GPU 2J 8� *K
fFRD
B__inference_cvae_3_layer_call_and_return_conditional_losses_1419642
StatefulPartitionedCall�
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*'
_output_shapes
:���������2

Identity"
identityIdentity:output:0*n
_input_shapes]
[:���������::::::::::::::::::22
StatefulPartitionedCallStatefulPartitionedCall:O K
'
_output_shapes
:���������
 
_user_specified_nameinputs
�
s
F__inference_sampling_3_layer_call_and_return_conditional_losses_141315

inputs
inputs_1
identity�D
ShapeShapeinputs*
T0*
_output_shapes
:2
Shapet
strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB: 2
strided_slice/stackx
strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice/stack_1x
strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice/stack_2�
strided_sliceStridedSliceShape:output:0strided_slice/stack:output:0strided_slice/stack_1:output:0strided_slice/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2
strided_sliceH
Shape_1Shapeinputs*
T0*
_output_shapes
:2	
Shape_1x
strided_slice_1/stackConst*
_output_shapes
:*
dtype0*
valueB:2
strided_slice_1/stack|
strided_slice_1/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice_1/stack_1|
strided_slice_1/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice_1/stack_2�
strided_slice_1StridedSliceShape_1:output:0strided_slice_1/stack:output:0 strided_slice_1/stack_1:output:0 strided_slice_1/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2
strided_slice_1�
random_normal/shapePackstrided_slice:output:0strided_slice_1:output:0*
N*
T0*
_output_shapes
:2
random_normal/shapem
random_normal/meanConst*
_output_shapes
: *
dtype0*
valueB
 *    2
random_normal/meanq
random_normal/stddevConst*
_output_shapes
: *
dtype0*
valueB
 *  �?2
random_normal/stddev�
"random_normal/RandomStandardNormalRandomStandardNormalrandom_normal/shape:output:0*
T0*0
_output_shapes
:������������������*
dtype0*
seed���)*
seed2���2$
"random_normal/RandomStandardNormal�
random_normal/mulMul+random_normal/RandomStandardNormal:output:0random_normal/stddev:output:0*
T0*0
_output_shapes
:������������������2
random_normal/mul�
random_normalAddrandom_normal/mul:z:0random_normal/mean:output:0*
T0*0
_output_shapes
:������������������2
random_normalS
mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *   ?2
mul/x]
mulMulmul/x:output:0inputs_1*
T0*'
_output_shapes
:���������2
mulL
ExpExpmul:z:0*
T0*'
_output_shapes
:���������2
Expc
mul_1MulExp:y:0random_normal:z:0*
T0*'
_output_shapes
:���������2
mul_1X
addAddV2inputs	mul_1:z:0*
T0*'
_output_shapes
:���������2
add[
IdentityIdentityadd:z:0*
T0*'
_output_shapes
:���������2

Identity"
identityIdentity:output:0*9
_input_shapes(
&:���������:���������:O K
'
_output_shapes
:���������
 
_user_specified_nameinputs:OK
'
_output_shapes
:���������
 
_user_specified_nameinputs
�
�
(__inference_encoder_layer_call_fn_141421
encoder_input
unknown
	unknown_0
	unknown_1
	unknown_2
	unknown_3
	unknown_4
	unknown_5
	unknown_6
	unknown_7
	unknown_8
identity

identity_1

identity_2��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallencoder_inputunknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6	unknown_7	unknown_8*
Tin
2*
Tout
2*
_collective_manager_ids
 *M
_output_shapes;
9:���������:���������:���������*,
_read_only_resource_inputs

	
*-
config_proto

CPU

GPU 2J 8� *L
fGRE
C__inference_encoder_layer_call_and_return_conditional_losses_1413942
StatefulPartitionedCall�
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*'
_output_shapes
:���������2

Identity�

Identity_1Identity StatefulPartitionedCall:output:1^StatefulPartitionedCall*
T0*'
_output_shapes
:���������2

Identity_1�

Identity_2Identity StatefulPartitionedCall:output:2^StatefulPartitionedCall*
T0*'
_output_shapes
:���������2

Identity_2"
identityIdentity:output:0"!

identity_1Identity_1:output:0"!

identity_2Identity_2:output:0*N
_input_shapes=
;:���������::::::::::22
StatefulPartitionedCallStatefulPartitionedCall:V R
'
_output_shapes
:���������
'
_user_specified_nameencoder_input
�
�
(__inference_decoder_layer_call_fn_141665
latent_variable
unknown
	unknown_0
	unknown_1
	unknown_2
	unknown_3
	unknown_4
	unknown_5
	unknown_6
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCalllatent_variableunknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6*
Tin
2	*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������**
_read_only_resource_inputs

*-
config_proto

CPU

GPU 2J 8� *L
fGRE
C__inference_decoder_layer_call_and_return_conditional_losses_1416462
StatefulPartitionedCall�
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*'
_output_shapes
:���������2

Identity"
identityIdentity:output:0*F
_input_shapes5
3:���������::::::::22
StatefulPartitionedCallStatefulPartitionedCall:X T
'
_output_shapes
:���������
)
_user_specified_namelatent_variable
�)
�
C__inference_decoder_layer_call_and_return_conditional_losses_142869

inputs-
)decoder_l2_matmul_readvariableop_resource.
*decoder_l2_biasadd_readvariableop_resource-
)decoder_l3_matmul_readvariableop_resource.
*decoder_l3_biasadd_readvariableop_resource-
)decoder_l4_matmul_readvariableop_resource.
*decoder_l4_biasadd_readvariableop_resource/
+output_layer_matmul_readvariableop_resource0
,output_layer_biasadd_readvariableop_resource
identity��!decoder_l2/BiasAdd/ReadVariableOp� decoder_l2/MatMul/ReadVariableOp�!decoder_l3/BiasAdd/ReadVariableOp� decoder_l3/MatMul/ReadVariableOp�!decoder_l4/BiasAdd/ReadVariableOp� decoder_l4/MatMul/ReadVariableOp�#output_layer/BiasAdd/ReadVariableOp�"output_layer/MatMul/ReadVariableOp�
 decoder_l2/MatMul/ReadVariableOpReadVariableOp)decoder_l2_matmul_readvariableop_resource*
_output_shapes

:*
dtype02"
 decoder_l2/MatMul/ReadVariableOp�
decoder_l2/MatMulMatMulinputs(decoder_l2/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2
decoder_l2/MatMul�
!decoder_l2/BiasAdd/ReadVariableOpReadVariableOp*decoder_l2_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02#
!decoder_l2/BiasAdd/ReadVariableOp�
decoder_l2/BiasAddBiasAdddecoder_l2/MatMul:product:0)decoder_l2/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2
decoder_l2/BiasAddy
decoder_l2/ReluReludecoder_l2/BiasAdd:output:0*
T0*'
_output_shapes
:���������2
decoder_l2/Relu�
 decoder_l3/MatMul/ReadVariableOpReadVariableOp)decoder_l3_matmul_readvariableop_resource*
_output_shapes

:*
dtype02"
 decoder_l3/MatMul/ReadVariableOp�
decoder_l3/MatMulMatMuldecoder_l2/Relu:activations:0(decoder_l3/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2
decoder_l3/MatMul�
!decoder_l3/BiasAdd/ReadVariableOpReadVariableOp*decoder_l3_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02#
!decoder_l3/BiasAdd/ReadVariableOp�
decoder_l3/BiasAddBiasAdddecoder_l3/MatMul:product:0)decoder_l3/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2
decoder_l3/BiasAddy
decoder_l3/ReluReludecoder_l3/BiasAdd:output:0*
T0*'
_output_shapes
:���������2
decoder_l3/Relu�
 decoder_l4/MatMul/ReadVariableOpReadVariableOp)decoder_l4_matmul_readvariableop_resource*
_output_shapes

:
*
dtype02"
 decoder_l4/MatMul/ReadVariableOp�
decoder_l4/MatMulMatMuldecoder_l3/Relu:activations:0(decoder_l4/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������
2
decoder_l4/MatMul�
!decoder_l4/BiasAdd/ReadVariableOpReadVariableOp*decoder_l4_biasadd_readvariableop_resource*
_output_shapes
:
*
dtype02#
!decoder_l4/BiasAdd/ReadVariableOp�
decoder_l4/BiasAddBiasAdddecoder_l4/MatMul:product:0)decoder_l4/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������
2
decoder_l4/BiasAddy
decoder_l4/ReluReludecoder_l4/BiasAdd:output:0*
T0*'
_output_shapes
:���������
2
decoder_l4/Relu�
"output_layer/MatMul/ReadVariableOpReadVariableOp+output_layer_matmul_readvariableop_resource*
_output_shapes

:
*
dtype02$
"output_layer/MatMul/ReadVariableOp�
output_layer/MatMulMatMuldecoder_l4/Relu:activations:0*output_layer/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2
output_layer/MatMul�
#output_layer/BiasAdd/ReadVariableOpReadVariableOp,output_layer_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02%
#output_layer/BiasAdd/ReadVariableOp�
output_layer/BiasAddBiasAddoutput_layer/MatMul:product:0+output_layer/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2
output_layer/BiasAdd
output_layer/ReluReluoutput_layer/BiasAdd:output:0*
T0*'
_output_shapes
:���������2
output_layer/Relu�
IdentityIdentityoutput_layer/Relu:activations:0"^decoder_l2/BiasAdd/ReadVariableOp!^decoder_l2/MatMul/ReadVariableOp"^decoder_l3/BiasAdd/ReadVariableOp!^decoder_l3/MatMul/ReadVariableOp"^decoder_l4/BiasAdd/ReadVariableOp!^decoder_l4/MatMul/ReadVariableOp$^output_layer/BiasAdd/ReadVariableOp#^output_layer/MatMul/ReadVariableOp*
T0*'
_output_shapes
:���������2

Identity"
identityIdentity:output:0*F
_input_shapes5
3:���������::::::::2F
!decoder_l2/BiasAdd/ReadVariableOp!decoder_l2/BiasAdd/ReadVariableOp2D
 decoder_l2/MatMul/ReadVariableOp decoder_l2/MatMul/ReadVariableOp2F
!decoder_l3/BiasAdd/ReadVariableOp!decoder_l3/BiasAdd/ReadVariableOp2D
 decoder_l3/MatMul/ReadVariableOp decoder_l3/MatMul/ReadVariableOp2F
!decoder_l4/BiasAdd/ReadVariableOp!decoder_l4/BiasAdd/ReadVariableOp2D
 decoder_l4/MatMul/ReadVariableOp decoder_l4/MatMul/ReadVariableOp2J
#output_layer/BiasAdd/ReadVariableOp#output_layer/BiasAdd/ReadVariableOp2H
"output_layer/MatMul/ReadVariableOp"output_layer/MatMul/ReadVariableOp:O K
'
_output_shapes
:���������
 
_user_specified_nameinputs
�	
�
F__inference_decoder_l2_layer_call_and_return_conditional_losses_143052

inputs"
matmul_readvariableop_resource#
biasadd_readvariableop_resource
identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOp�
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:*
dtype02
MatMul/ReadVariableOps
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2
MatMul�
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:*
dtype02
BiasAdd/ReadVariableOp�
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2	
BiasAddX
ReluReluBiasAdd:output:0*
T0*'
_output_shapes
:���������2
Relu�
IdentityIdentityRelu:activations:0^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*
T0*'
_output_shapes
:���������2

Identity"
identityIdentity:output:0*.
_input_shapes
:���������::20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:O K
'
_output_shapes
:���������
 
_user_specified_nameinputs
�	
�
F__inference_decoder_l4_layer_call_and_return_conditional_losses_141551

inputs"
matmul_readvariableop_resource#
biasadd_readvariableop_resource
identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOp�
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:
*
dtype02
MatMul/ReadVariableOps
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������
2
MatMul�
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:
*
dtype02
BiasAdd/ReadVariableOp�
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������
2	
BiasAddX
ReluReluBiasAdd:output:0*
T0*'
_output_shapes
:���������
2
Relu�
IdentityIdentityRelu:activations:0^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*
T0*'
_output_shapes
:���������
2

Identity"
identityIdentity:output:0*.
_input_shapes
:���������::20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:O K
'
_output_shapes
:���������
 
_user_specified_nameinputs
��
�
B__inference_cvae_3_layer_call_and_return_conditional_losses_142279
input_15
1encoder_encoder_l2_matmul_readvariableop_resource6
2encoder_encoder_l2_biasadd_readvariableop_resource5
1encoder_encoder_l3_matmul_readvariableop_resource6
2encoder_encoder_l3_biasadd_readvariableop_resource5
1encoder_encoder_l4_matmul_readvariableop_resource6
2encoder_encoder_l4_biasadd_readvariableop_resource1
-encoder_z_mean_matmul_readvariableop_resource2
.encoder_z_mean_biasadd_readvariableop_resource4
0encoder_z_log_var_matmul_readvariableop_resource5
1encoder_z_log_var_biasadd_readvariableop_resource5
1decoder_decoder_l2_matmul_readvariableop_resource6
2decoder_decoder_l2_biasadd_readvariableop_resource5
1decoder_decoder_l3_matmul_readvariableop_resource6
2decoder_decoder_l3_biasadd_readvariableop_resource5
1decoder_decoder_l4_matmul_readvariableop_resource6
2decoder_decoder_l4_biasadd_readvariableop_resource7
3decoder_output_layer_matmul_readvariableop_resource8
4decoder_output_layer_biasadd_readvariableop_resource
identity��)decoder/decoder_l2/BiasAdd/ReadVariableOp�(decoder/decoder_l2/MatMul/ReadVariableOp�)decoder/decoder_l3/BiasAdd/ReadVariableOp�(decoder/decoder_l3/MatMul/ReadVariableOp�)decoder/decoder_l4/BiasAdd/ReadVariableOp�(decoder/decoder_l4/MatMul/ReadVariableOp�+decoder/output_layer/BiasAdd/ReadVariableOp�*decoder/output_layer/MatMul/ReadVariableOp�)encoder/encoder_l2/BiasAdd/ReadVariableOp�(encoder/encoder_l2/MatMul/ReadVariableOp�)encoder/encoder_l3/BiasAdd/ReadVariableOp�(encoder/encoder_l3/MatMul/ReadVariableOp�)encoder/encoder_l4/BiasAdd/ReadVariableOp�(encoder/encoder_l4/MatMul/ReadVariableOp�(encoder/z_log_var/BiasAdd/ReadVariableOp�'encoder/z_log_var/MatMul/ReadVariableOp�%encoder/z_mean/BiasAdd/ReadVariableOp�$encoder/z_mean/MatMul/ReadVariableOp�
(encoder/encoder_l2/MatMul/ReadVariableOpReadVariableOp1encoder_encoder_l2_matmul_readvariableop_resource*
_output_shapes

:
*
dtype02*
(encoder/encoder_l2/MatMul/ReadVariableOp�
encoder/encoder_l2/MatMulMatMulinput_10encoder/encoder_l2/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������
2
encoder/encoder_l2/MatMul�
)encoder/encoder_l2/BiasAdd/ReadVariableOpReadVariableOp2encoder_encoder_l2_biasadd_readvariableop_resource*
_output_shapes
:
*
dtype02+
)encoder/encoder_l2/BiasAdd/ReadVariableOp�
encoder/encoder_l2/BiasAddBiasAdd#encoder/encoder_l2/MatMul:product:01encoder/encoder_l2/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������
2
encoder/encoder_l2/BiasAdd�
encoder/encoder_l2/ReluRelu#encoder/encoder_l2/BiasAdd:output:0*
T0*'
_output_shapes
:���������
2
encoder/encoder_l2/Relu�
(encoder/encoder_l3/MatMul/ReadVariableOpReadVariableOp1encoder_encoder_l3_matmul_readvariableop_resource*
_output_shapes

:
*
dtype02*
(encoder/encoder_l3/MatMul/ReadVariableOp�
encoder/encoder_l3/MatMulMatMul%encoder/encoder_l2/Relu:activations:00encoder/encoder_l3/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2
encoder/encoder_l3/MatMul�
)encoder/encoder_l3/BiasAdd/ReadVariableOpReadVariableOp2encoder_encoder_l3_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02+
)encoder/encoder_l3/BiasAdd/ReadVariableOp�
encoder/encoder_l3/BiasAddBiasAdd#encoder/encoder_l3/MatMul:product:01encoder/encoder_l3/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2
encoder/encoder_l3/BiasAdd�
encoder/encoder_l3/ReluRelu#encoder/encoder_l3/BiasAdd:output:0*
T0*'
_output_shapes
:���������2
encoder/encoder_l3/Relu�
(encoder/encoder_l4/MatMul/ReadVariableOpReadVariableOp1encoder_encoder_l4_matmul_readvariableop_resource*
_output_shapes

:*
dtype02*
(encoder/encoder_l4/MatMul/ReadVariableOp�
encoder/encoder_l4/MatMulMatMul%encoder/encoder_l3/Relu:activations:00encoder/encoder_l4/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2
encoder/encoder_l4/MatMul�
)encoder/encoder_l4/BiasAdd/ReadVariableOpReadVariableOp2encoder_encoder_l4_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02+
)encoder/encoder_l4/BiasAdd/ReadVariableOp�
encoder/encoder_l4/BiasAddBiasAdd#encoder/encoder_l4/MatMul:product:01encoder/encoder_l4/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2
encoder/encoder_l4/BiasAdd�
encoder/encoder_l4/ReluRelu#encoder/encoder_l4/BiasAdd:output:0*
T0*'
_output_shapes
:���������2
encoder/encoder_l4/Relu�
$encoder/z_mean/MatMul/ReadVariableOpReadVariableOp-encoder_z_mean_matmul_readvariableop_resource*
_output_shapes

:*
dtype02&
$encoder/z_mean/MatMul/ReadVariableOp�
encoder/z_mean/MatMulMatMul%encoder/encoder_l4/Relu:activations:0,encoder/z_mean/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2
encoder/z_mean/MatMul�
%encoder/z_mean/BiasAdd/ReadVariableOpReadVariableOp.encoder_z_mean_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02'
%encoder/z_mean/BiasAdd/ReadVariableOp�
encoder/z_mean/BiasAddBiasAddencoder/z_mean/MatMul:product:0-encoder/z_mean/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2
encoder/z_mean/BiasAdd�
'encoder/z_log_var/MatMul/ReadVariableOpReadVariableOp0encoder_z_log_var_matmul_readvariableop_resource*
_output_shapes

:*
dtype02)
'encoder/z_log_var/MatMul/ReadVariableOp�
encoder/z_log_var/MatMulMatMul%encoder/encoder_l4/Relu:activations:0/encoder/z_log_var/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2
encoder/z_log_var/MatMul�
(encoder/z_log_var/BiasAdd/ReadVariableOpReadVariableOp1encoder_z_log_var_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02*
(encoder/z_log_var/BiasAdd/ReadVariableOp�
encoder/z_log_var/BiasAddBiasAdd"encoder/z_log_var/MatMul:product:00encoder/z_log_var/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2
encoder/z_log_var/BiasAdd�
encoder/sampling_3/ShapeShapeencoder/z_mean/BiasAdd:output:0*
T0*
_output_shapes
:2
encoder/sampling_3/Shape�
&encoder/sampling_3/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB: 2(
&encoder/sampling_3/strided_slice/stack�
(encoder/sampling_3/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2*
(encoder/sampling_3/strided_slice/stack_1�
(encoder/sampling_3/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2*
(encoder/sampling_3/strided_slice/stack_2�
 encoder/sampling_3/strided_sliceStridedSlice!encoder/sampling_3/Shape:output:0/encoder/sampling_3/strided_slice/stack:output:01encoder/sampling_3/strided_slice/stack_1:output:01encoder/sampling_3/strided_slice/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2"
 encoder/sampling_3/strided_slice�
encoder/sampling_3/Shape_1Shapeencoder/z_mean/BiasAdd:output:0*
T0*
_output_shapes
:2
encoder/sampling_3/Shape_1�
(encoder/sampling_3/strided_slice_1/stackConst*
_output_shapes
:*
dtype0*
valueB:2*
(encoder/sampling_3/strided_slice_1/stack�
*encoder/sampling_3/strided_slice_1/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2,
*encoder/sampling_3/strided_slice_1/stack_1�
*encoder/sampling_3/strided_slice_1/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2,
*encoder/sampling_3/strided_slice_1/stack_2�
"encoder/sampling_3/strided_slice_1StridedSlice#encoder/sampling_3/Shape_1:output:01encoder/sampling_3/strided_slice_1/stack:output:03encoder/sampling_3/strided_slice_1/stack_1:output:03encoder/sampling_3/strided_slice_1/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2$
"encoder/sampling_3/strided_slice_1�
&encoder/sampling_3/random_normal/shapePack)encoder/sampling_3/strided_slice:output:0+encoder/sampling_3/strided_slice_1:output:0*
N*
T0*
_output_shapes
:2(
&encoder/sampling_3/random_normal/shape�
%encoder/sampling_3/random_normal/meanConst*
_output_shapes
: *
dtype0*
valueB
 *    2'
%encoder/sampling_3/random_normal/mean�
'encoder/sampling_3/random_normal/stddevConst*
_output_shapes
: *
dtype0*
valueB
 *  �?2)
'encoder/sampling_3/random_normal/stddev�
5encoder/sampling_3/random_normal/RandomStandardNormalRandomStandardNormal/encoder/sampling_3/random_normal/shape:output:0*
T0*0
_output_shapes
:������������������*
dtype0*
seed���)*
seed2���27
5encoder/sampling_3/random_normal/RandomStandardNormal�
$encoder/sampling_3/random_normal/mulMul>encoder/sampling_3/random_normal/RandomStandardNormal:output:00encoder/sampling_3/random_normal/stddev:output:0*
T0*0
_output_shapes
:������������������2&
$encoder/sampling_3/random_normal/mul�
 encoder/sampling_3/random_normalAdd(encoder/sampling_3/random_normal/mul:z:0.encoder/sampling_3/random_normal/mean:output:0*
T0*0
_output_shapes
:������������������2"
 encoder/sampling_3/random_normaly
encoder/sampling_3/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *   ?2
encoder/sampling_3/mul/x�
encoder/sampling_3/mulMul!encoder/sampling_3/mul/x:output:0"encoder/z_log_var/BiasAdd:output:0*
T0*'
_output_shapes
:���������2
encoder/sampling_3/mul�
encoder/sampling_3/ExpExpencoder/sampling_3/mul:z:0*
T0*'
_output_shapes
:���������2
encoder/sampling_3/Exp�
encoder/sampling_3/mul_1Mulencoder/sampling_3/Exp:y:0$encoder/sampling_3/random_normal:z:0*
T0*'
_output_shapes
:���������2
encoder/sampling_3/mul_1�
encoder/sampling_3/addAddV2encoder/z_mean/BiasAdd:output:0encoder/sampling_3/mul_1:z:0*
T0*'
_output_shapes
:���������2
encoder/sampling_3/add{
strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB"    ����2
strided_slice/stack
strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB"        2
strided_slice/stack_1
strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      2
strided_slice/stack_2�
strided_sliceStridedSliceinput_1strided_slice/stack:output:0strided_slice/stack_1:output:0strided_slice/stack_2:output:0*
Index0*
T0*'
_output_shapes
:���������*

begin_mask*
end_mask2
strided_slicet
concatenate/concat/axisConst*
_output_shapes
: *
dtype0*
value	B :2
concatenate/concat/axis�
concatenate/concatConcatV2encoder/sampling_3/add:z:0strided_slice:output:0 concatenate/concat/axis:output:0*
N*
T0*'
_output_shapes
:���������2
concatenate/concat�
(decoder/decoder_l2/MatMul/ReadVariableOpReadVariableOp1decoder_decoder_l2_matmul_readvariableop_resource*
_output_shapes

:*
dtype02*
(decoder/decoder_l2/MatMul/ReadVariableOp�
decoder/decoder_l2/MatMulMatMulconcatenate/concat:output:00decoder/decoder_l2/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2
decoder/decoder_l2/MatMul�
)decoder/decoder_l2/BiasAdd/ReadVariableOpReadVariableOp2decoder_decoder_l2_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02+
)decoder/decoder_l2/BiasAdd/ReadVariableOp�
decoder/decoder_l2/BiasAddBiasAdd#decoder/decoder_l2/MatMul:product:01decoder/decoder_l2/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2
decoder/decoder_l2/BiasAdd�
decoder/decoder_l2/ReluRelu#decoder/decoder_l2/BiasAdd:output:0*
T0*'
_output_shapes
:���������2
decoder/decoder_l2/Relu�
(decoder/decoder_l3/MatMul/ReadVariableOpReadVariableOp1decoder_decoder_l3_matmul_readvariableop_resource*
_output_shapes

:*
dtype02*
(decoder/decoder_l3/MatMul/ReadVariableOp�
decoder/decoder_l3/MatMulMatMul%decoder/decoder_l2/Relu:activations:00decoder/decoder_l3/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2
decoder/decoder_l3/MatMul�
)decoder/decoder_l3/BiasAdd/ReadVariableOpReadVariableOp2decoder_decoder_l3_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02+
)decoder/decoder_l3/BiasAdd/ReadVariableOp�
decoder/decoder_l3/BiasAddBiasAdd#decoder/decoder_l3/MatMul:product:01decoder/decoder_l3/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2
decoder/decoder_l3/BiasAdd�
decoder/decoder_l3/ReluRelu#decoder/decoder_l3/BiasAdd:output:0*
T0*'
_output_shapes
:���������2
decoder/decoder_l3/Relu�
(decoder/decoder_l4/MatMul/ReadVariableOpReadVariableOp1decoder_decoder_l4_matmul_readvariableop_resource*
_output_shapes

:
*
dtype02*
(decoder/decoder_l4/MatMul/ReadVariableOp�
decoder/decoder_l4/MatMulMatMul%decoder/decoder_l3/Relu:activations:00decoder/decoder_l4/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������
2
decoder/decoder_l4/MatMul�
)decoder/decoder_l4/BiasAdd/ReadVariableOpReadVariableOp2decoder_decoder_l4_biasadd_readvariableop_resource*
_output_shapes
:
*
dtype02+
)decoder/decoder_l4/BiasAdd/ReadVariableOp�
decoder/decoder_l4/BiasAddBiasAdd#decoder/decoder_l4/MatMul:product:01decoder/decoder_l4/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������
2
decoder/decoder_l4/BiasAdd�
decoder/decoder_l4/ReluRelu#decoder/decoder_l4/BiasAdd:output:0*
T0*'
_output_shapes
:���������
2
decoder/decoder_l4/Relu�
*decoder/output_layer/MatMul/ReadVariableOpReadVariableOp3decoder_output_layer_matmul_readvariableop_resource*
_output_shapes

:
*
dtype02,
*decoder/output_layer/MatMul/ReadVariableOp�
decoder/output_layer/MatMulMatMul%decoder/decoder_l4/Relu:activations:02decoder/output_layer/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2
decoder/output_layer/MatMul�
+decoder/output_layer/BiasAdd/ReadVariableOpReadVariableOp4decoder_output_layer_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02-
+decoder/output_layer/BiasAdd/ReadVariableOp�
decoder/output_layer/BiasAddBiasAdd%decoder/output_layer/MatMul:product:03decoder/output_layer/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2
decoder/output_layer/BiasAdd�
decoder/output_layer/ReluRelu%decoder/output_layer/BiasAdd:output:0*
T0*'
_output_shapes
:���������2
decoder/output_layer/Relu�
IdentityIdentity'decoder/output_layer/Relu:activations:0*^decoder/decoder_l2/BiasAdd/ReadVariableOp)^decoder/decoder_l2/MatMul/ReadVariableOp*^decoder/decoder_l3/BiasAdd/ReadVariableOp)^decoder/decoder_l3/MatMul/ReadVariableOp*^decoder/decoder_l4/BiasAdd/ReadVariableOp)^decoder/decoder_l4/MatMul/ReadVariableOp,^decoder/output_layer/BiasAdd/ReadVariableOp+^decoder/output_layer/MatMul/ReadVariableOp*^encoder/encoder_l2/BiasAdd/ReadVariableOp)^encoder/encoder_l2/MatMul/ReadVariableOp*^encoder/encoder_l3/BiasAdd/ReadVariableOp)^encoder/encoder_l3/MatMul/ReadVariableOp*^encoder/encoder_l4/BiasAdd/ReadVariableOp)^encoder/encoder_l4/MatMul/ReadVariableOp)^encoder/z_log_var/BiasAdd/ReadVariableOp(^encoder/z_log_var/MatMul/ReadVariableOp&^encoder/z_mean/BiasAdd/ReadVariableOp%^encoder/z_mean/MatMul/ReadVariableOp*
T0*'
_output_shapes
:���������2

Identity"
identityIdentity:output:0*n
_input_shapes]
[:���������::::::::::::::::::2V
)decoder/decoder_l2/BiasAdd/ReadVariableOp)decoder/decoder_l2/BiasAdd/ReadVariableOp2T
(decoder/decoder_l2/MatMul/ReadVariableOp(decoder/decoder_l2/MatMul/ReadVariableOp2V
)decoder/decoder_l3/BiasAdd/ReadVariableOp)decoder/decoder_l3/BiasAdd/ReadVariableOp2T
(decoder/decoder_l3/MatMul/ReadVariableOp(decoder/decoder_l3/MatMul/ReadVariableOp2V
)decoder/decoder_l4/BiasAdd/ReadVariableOp)decoder/decoder_l4/BiasAdd/ReadVariableOp2T
(decoder/decoder_l4/MatMul/ReadVariableOp(decoder/decoder_l4/MatMul/ReadVariableOp2Z
+decoder/output_layer/BiasAdd/ReadVariableOp+decoder/output_layer/BiasAdd/ReadVariableOp2X
*decoder/output_layer/MatMul/ReadVariableOp*decoder/output_layer/MatMul/ReadVariableOp2V
)encoder/encoder_l2/BiasAdd/ReadVariableOp)encoder/encoder_l2/BiasAdd/ReadVariableOp2T
(encoder/encoder_l2/MatMul/ReadVariableOp(encoder/encoder_l2/MatMul/ReadVariableOp2V
)encoder/encoder_l3/BiasAdd/ReadVariableOp)encoder/encoder_l3/BiasAdd/ReadVariableOp2T
(encoder/encoder_l3/MatMul/ReadVariableOp(encoder/encoder_l3/MatMul/ReadVariableOp2V
)encoder/encoder_l4/BiasAdd/ReadVariableOp)encoder/encoder_l4/BiasAdd/ReadVariableOp2T
(encoder/encoder_l4/MatMul/ReadVariableOp(encoder/encoder_l4/MatMul/ReadVariableOp2T
(encoder/z_log_var/BiasAdd/ReadVariableOp(encoder/z_log_var/BiasAdd/ReadVariableOp2R
'encoder/z_log_var/MatMul/ReadVariableOp'encoder/z_log_var/MatMul/ReadVariableOp2N
%encoder/z_mean/BiasAdd/ReadVariableOp%encoder/z_mean/BiasAdd/ReadVariableOp2L
$encoder/z_mean/MatMul/ReadVariableOp$encoder/z_mean/MatMul/ReadVariableOp:P L
'
_output_shapes
:���������
!
_user_specified_name	input_1
�	
�
F__inference_encoder_l3_layer_call_and_return_conditional_losses_142942

inputs"
matmul_readvariableop_resource#
biasadd_readvariableop_resource
identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOp�
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:
*
dtype02
MatMul/ReadVariableOps
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2
MatMul�
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:*
dtype02
BiasAdd/ReadVariableOp�
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2	
BiasAddX
ReluReluBiasAdd:output:0*
T0*'
_output_shapes
:���������2
Relu�
IdentityIdentityRelu:activations:0^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*
T0*'
_output_shapes
:���������2

Identity"
identityIdentity:output:0*.
_input_shapes
:���������
::20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:O K
'
_output_shapes
:���������

 
_user_specified_nameinputs
�

�
'__inference_cvae_3_layer_call_fn_142320
input_1
unknown
	unknown_0
	unknown_1
	unknown_2
	unknown_3
	unknown_4
	unknown_5
	unknown_6
	unknown_7
	unknown_8
	unknown_9

unknown_10

unknown_11

unknown_12

unknown_13

unknown_14

unknown_15

unknown_16
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinput_1unknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6	unknown_7	unknown_8	unknown_9
unknown_10
unknown_11
unknown_12
unknown_13
unknown_14
unknown_15
unknown_16*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*4
_read_only_resource_inputs
	
*-
config_proto

CPU

GPU 2J 8� *K
fFRD
B__inference_cvae_3_layer_call_and_return_conditional_losses_1419642
StatefulPartitionedCall�
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*'
_output_shapes
:���������2

Identity"
identityIdentity:output:0*n
_input_shapes]
[:���������::::::::::::::::::22
StatefulPartitionedCallStatefulPartitionedCall:P L
'
_output_shapes
:���������
!
_user_specified_name	input_1
�)
�
C__inference_decoder_layer_call_and_return_conditional_losses_142837

inputs-
)decoder_l2_matmul_readvariableop_resource.
*decoder_l2_biasadd_readvariableop_resource-
)decoder_l3_matmul_readvariableop_resource.
*decoder_l3_biasadd_readvariableop_resource-
)decoder_l4_matmul_readvariableop_resource.
*decoder_l4_biasadd_readvariableop_resource/
+output_layer_matmul_readvariableop_resource0
,output_layer_biasadd_readvariableop_resource
identity��!decoder_l2/BiasAdd/ReadVariableOp� decoder_l2/MatMul/ReadVariableOp�!decoder_l3/BiasAdd/ReadVariableOp� decoder_l3/MatMul/ReadVariableOp�!decoder_l4/BiasAdd/ReadVariableOp� decoder_l4/MatMul/ReadVariableOp�#output_layer/BiasAdd/ReadVariableOp�"output_layer/MatMul/ReadVariableOp�
 decoder_l2/MatMul/ReadVariableOpReadVariableOp)decoder_l2_matmul_readvariableop_resource*
_output_shapes

:*
dtype02"
 decoder_l2/MatMul/ReadVariableOp�
decoder_l2/MatMulMatMulinputs(decoder_l2/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2
decoder_l2/MatMul�
!decoder_l2/BiasAdd/ReadVariableOpReadVariableOp*decoder_l2_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02#
!decoder_l2/BiasAdd/ReadVariableOp�
decoder_l2/BiasAddBiasAdddecoder_l2/MatMul:product:0)decoder_l2/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2
decoder_l2/BiasAddy
decoder_l2/ReluReludecoder_l2/BiasAdd:output:0*
T0*'
_output_shapes
:���������2
decoder_l2/Relu�
 decoder_l3/MatMul/ReadVariableOpReadVariableOp)decoder_l3_matmul_readvariableop_resource*
_output_shapes

:*
dtype02"
 decoder_l3/MatMul/ReadVariableOp�
decoder_l3/MatMulMatMuldecoder_l2/Relu:activations:0(decoder_l3/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2
decoder_l3/MatMul�
!decoder_l3/BiasAdd/ReadVariableOpReadVariableOp*decoder_l3_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02#
!decoder_l3/BiasAdd/ReadVariableOp�
decoder_l3/BiasAddBiasAdddecoder_l3/MatMul:product:0)decoder_l3/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2
decoder_l3/BiasAddy
decoder_l3/ReluReludecoder_l3/BiasAdd:output:0*
T0*'
_output_shapes
:���������2
decoder_l3/Relu�
 decoder_l4/MatMul/ReadVariableOpReadVariableOp)decoder_l4_matmul_readvariableop_resource*
_output_shapes

:
*
dtype02"
 decoder_l4/MatMul/ReadVariableOp�
decoder_l4/MatMulMatMuldecoder_l3/Relu:activations:0(decoder_l4/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������
2
decoder_l4/MatMul�
!decoder_l4/BiasAdd/ReadVariableOpReadVariableOp*decoder_l4_biasadd_readvariableop_resource*
_output_shapes
:
*
dtype02#
!decoder_l4/BiasAdd/ReadVariableOp�
decoder_l4/BiasAddBiasAdddecoder_l4/MatMul:product:0)decoder_l4/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������
2
decoder_l4/BiasAddy
decoder_l4/ReluReludecoder_l4/BiasAdd:output:0*
T0*'
_output_shapes
:���������
2
decoder_l4/Relu�
"output_layer/MatMul/ReadVariableOpReadVariableOp+output_layer_matmul_readvariableop_resource*
_output_shapes

:
*
dtype02$
"output_layer/MatMul/ReadVariableOp�
output_layer/MatMulMatMuldecoder_l4/Relu:activations:0*output_layer/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2
output_layer/MatMul�
#output_layer/BiasAdd/ReadVariableOpReadVariableOp,output_layer_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02%
#output_layer/BiasAdd/ReadVariableOp�
output_layer/BiasAddBiasAddoutput_layer/MatMul:product:0+output_layer/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2
output_layer/BiasAdd
output_layer/ReluReluoutput_layer/BiasAdd:output:0*
T0*'
_output_shapes
:���������2
output_layer/Relu�
IdentityIdentityoutput_layer/Relu:activations:0"^decoder_l2/BiasAdd/ReadVariableOp!^decoder_l2/MatMul/ReadVariableOp"^decoder_l3/BiasAdd/ReadVariableOp!^decoder_l3/MatMul/ReadVariableOp"^decoder_l4/BiasAdd/ReadVariableOp!^decoder_l4/MatMul/ReadVariableOp$^output_layer/BiasAdd/ReadVariableOp#^output_layer/MatMul/ReadVariableOp*
T0*'
_output_shapes
:���������2

Identity"
identityIdentity:output:0*F
_input_shapes5
3:���������::::::::2F
!decoder_l2/BiasAdd/ReadVariableOp!decoder_l2/BiasAdd/ReadVariableOp2D
 decoder_l2/MatMul/ReadVariableOp decoder_l2/MatMul/ReadVariableOp2F
!decoder_l3/BiasAdd/ReadVariableOp!decoder_l3/BiasAdd/ReadVariableOp2D
 decoder_l3/MatMul/ReadVariableOp decoder_l3/MatMul/ReadVariableOp2F
!decoder_l4/BiasAdd/ReadVariableOp!decoder_l4/BiasAdd/ReadVariableOp2D
 decoder_l4/MatMul/ReadVariableOp decoder_l4/MatMul/ReadVariableOp2J
#output_layer/BiasAdd/ReadVariableOp#output_layer/BiasAdd/ReadVariableOp2H
"output_layer/MatMul/ReadVariableOp"output_layer/MatMul/ReadVariableOp:O K
'
_output_shapes
:���������
 
_user_specified_nameinputs
�
�
-__inference_output_layer_layer_call_fn_143121

inputs
unknown
	unknown_0
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *Q
fLRJ
H__inference_output_layer_layer_call_and_return_conditional_losses_1415782
StatefulPartitionedCall�
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*'
_output_shapes
:���������2

Identity"
identityIdentity:output:0*.
_input_shapes
:���������
::22
StatefulPartitionedCallStatefulPartitionedCall:O K
'
_output_shapes
:���������

 
_user_specified_nameinputs
�	
�
F__inference_encoder_l4_layer_call_and_return_conditional_losses_141221

inputs"
matmul_readvariableop_resource#
biasadd_readvariableop_resource
identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOp�
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:*
dtype02
MatMul/ReadVariableOps
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2
MatMul�
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:*
dtype02
BiasAdd/ReadVariableOp�
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2	
BiasAddX
ReluReluBiasAdd:output:0*
T0*'
_output_shapes
:���������2
Relu�
IdentityIdentityRelu:activations:0^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*
T0*'
_output_shapes
:���������2

Identity"
identityIdentity:output:0*.
_input_shapes
:���������::20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:O K
'
_output_shapes
:���������
 
_user_specified_nameinputs
�	
�
F__inference_encoder_l2_layer_call_and_return_conditional_losses_142922

inputs"
matmul_readvariableop_resource#
biasadd_readvariableop_resource
identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOp�
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:
*
dtype02
MatMul/ReadVariableOps
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������
2
MatMul�
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:
*
dtype02
BiasAdd/ReadVariableOp�
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������
2	
BiasAddX
ReluReluBiasAdd:output:0*
T0*'
_output_shapes
:���������
2
Relu�
IdentityIdentityRelu:activations:0^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*
T0*'
_output_shapes
:���������
2

Identity"
identityIdentity:output:0*.
_input_shapes
:���������::20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:O K
'
_output_shapes
:���������
 
_user_specified_nameinputs
��
�
B__inference_cvae_3_layer_call_and_return_conditional_losses_142453

inputs5
1encoder_encoder_l2_matmul_readvariableop_resource6
2encoder_encoder_l2_biasadd_readvariableop_resource5
1encoder_encoder_l3_matmul_readvariableop_resource6
2encoder_encoder_l3_biasadd_readvariableop_resource5
1encoder_encoder_l4_matmul_readvariableop_resource6
2encoder_encoder_l4_biasadd_readvariableop_resource1
-encoder_z_mean_matmul_readvariableop_resource2
.encoder_z_mean_biasadd_readvariableop_resource4
0encoder_z_log_var_matmul_readvariableop_resource5
1encoder_z_log_var_biasadd_readvariableop_resource5
1decoder_decoder_l2_matmul_readvariableop_resource6
2decoder_decoder_l2_biasadd_readvariableop_resource5
1decoder_decoder_l3_matmul_readvariableop_resource6
2decoder_decoder_l3_biasadd_readvariableop_resource5
1decoder_decoder_l4_matmul_readvariableop_resource6
2decoder_decoder_l4_biasadd_readvariableop_resource7
3decoder_output_layer_matmul_readvariableop_resource8
4decoder_output_layer_biasadd_readvariableop_resource
identity��)decoder/decoder_l2/BiasAdd/ReadVariableOp�(decoder/decoder_l2/MatMul/ReadVariableOp�)decoder/decoder_l3/BiasAdd/ReadVariableOp�(decoder/decoder_l3/MatMul/ReadVariableOp�)decoder/decoder_l4/BiasAdd/ReadVariableOp�(decoder/decoder_l4/MatMul/ReadVariableOp�+decoder/output_layer/BiasAdd/ReadVariableOp�*decoder/output_layer/MatMul/ReadVariableOp�)encoder/encoder_l2/BiasAdd/ReadVariableOp�(encoder/encoder_l2/MatMul/ReadVariableOp�)encoder/encoder_l3/BiasAdd/ReadVariableOp�(encoder/encoder_l3/MatMul/ReadVariableOp�)encoder/encoder_l4/BiasAdd/ReadVariableOp�(encoder/encoder_l4/MatMul/ReadVariableOp�(encoder/z_log_var/BiasAdd/ReadVariableOp�'encoder/z_log_var/MatMul/ReadVariableOp�%encoder/z_mean/BiasAdd/ReadVariableOp�$encoder/z_mean/MatMul/ReadVariableOp�
(encoder/encoder_l2/MatMul/ReadVariableOpReadVariableOp1encoder_encoder_l2_matmul_readvariableop_resource*
_output_shapes

:
*
dtype02*
(encoder/encoder_l2/MatMul/ReadVariableOp�
encoder/encoder_l2/MatMulMatMulinputs0encoder/encoder_l2/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������
2
encoder/encoder_l2/MatMul�
)encoder/encoder_l2/BiasAdd/ReadVariableOpReadVariableOp2encoder_encoder_l2_biasadd_readvariableop_resource*
_output_shapes
:
*
dtype02+
)encoder/encoder_l2/BiasAdd/ReadVariableOp�
encoder/encoder_l2/BiasAddBiasAdd#encoder/encoder_l2/MatMul:product:01encoder/encoder_l2/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������
2
encoder/encoder_l2/BiasAdd�
encoder/encoder_l2/ReluRelu#encoder/encoder_l2/BiasAdd:output:0*
T0*'
_output_shapes
:���������
2
encoder/encoder_l2/Relu�
(encoder/encoder_l3/MatMul/ReadVariableOpReadVariableOp1encoder_encoder_l3_matmul_readvariableop_resource*
_output_shapes

:
*
dtype02*
(encoder/encoder_l3/MatMul/ReadVariableOp�
encoder/encoder_l3/MatMulMatMul%encoder/encoder_l2/Relu:activations:00encoder/encoder_l3/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2
encoder/encoder_l3/MatMul�
)encoder/encoder_l3/BiasAdd/ReadVariableOpReadVariableOp2encoder_encoder_l3_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02+
)encoder/encoder_l3/BiasAdd/ReadVariableOp�
encoder/encoder_l3/BiasAddBiasAdd#encoder/encoder_l3/MatMul:product:01encoder/encoder_l3/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2
encoder/encoder_l3/BiasAdd�
encoder/encoder_l3/ReluRelu#encoder/encoder_l3/BiasAdd:output:0*
T0*'
_output_shapes
:���������2
encoder/encoder_l3/Relu�
(encoder/encoder_l4/MatMul/ReadVariableOpReadVariableOp1encoder_encoder_l4_matmul_readvariableop_resource*
_output_shapes

:*
dtype02*
(encoder/encoder_l4/MatMul/ReadVariableOp�
encoder/encoder_l4/MatMulMatMul%encoder/encoder_l3/Relu:activations:00encoder/encoder_l4/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2
encoder/encoder_l4/MatMul�
)encoder/encoder_l4/BiasAdd/ReadVariableOpReadVariableOp2encoder_encoder_l4_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02+
)encoder/encoder_l4/BiasAdd/ReadVariableOp�
encoder/encoder_l4/BiasAddBiasAdd#encoder/encoder_l4/MatMul:product:01encoder/encoder_l4/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2
encoder/encoder_l4/BiasAdd�
encoder/encoder_l4/ReluRelu#encoder/encoder_l4/BiasAdd:output:0*
T0*'
_output_shapes
:���������2
encoder/encoder_l4/Relu�
$encoder/z_mean/MatMul/ReadVariableOpReadVariableOp-encoder_z_mean_matmul_readvariableop_resource*
_output_shapes

:*
dtype02&
$encoder/z_mean/MatMul/ReadVariableOp�
encoder/z_mean/MatMulMatMul%encoder/encoder_l4/Relu:activations:0,encoder/z_mean/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2
encoder/z_mean/MatMul�
%encoder/z_mean/BiasAdd/ReadVariableOpReadVariableOp.encoder_z_mean_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02'
%encoder/z_mean/BiasAdd/ReadVariableOp�
encoder/z_mean/BiasAddBiasAddencoder/z_mean/MatMul:product:0-encoder/z_mean/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2
encoder/z_mean/BiasAdd�
'encoder/z_log_var/MatMul/ReadVariableOpReadVariableOp0encoder_z_log_var_matmul_readvariableop_resource*
_output_shapes

:*
dtype02)
'encoder/z_log_var/MatMul/ReadVariableOp�
encoder/z_log_var/MatMulMatMul%encoder/encoder_l4/Relu:activations:0/encoder/z_log_var/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2
encoder/z_log_var/MatMul�
(encoder/z_log_var/BiasAdd/ReadVariableOpReadVariableOp1encoder_z_log_var_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02*
(encoder/z_log_var/BiasAdd/ReadVariableOp�
encoder/z_log_var/BiasAddBiasAdd"encoder/z_log_var/MatMul:product:00encoder/z_log_var/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2
encoder/z_log_var/BiasAdd�
encoder/sampling_3/ShapeShapeencoder/z_mean/BiasAdd:output:0*
T0*
_output_shapes
:2
encoder/sampling_3/Shape�
&encoder/sampling_3/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB: 2(
&encoder/sampling_3/strided_slice/stack�
(encoder/sampling_3/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2*
(encoder/sampling_3/strided_slice/stack_1�
(encoder/sampling_3/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2*
(encoder/sampling_3/strided_slice/stack_2�
 encoder/sampling_3/strided_sliceStridedSlice!encoder/sampling_3/Shape:output:0/encoder/sampling_3/strided_slice/stack:output:01encoder/sampling_3/strided_slice/stack_1:output:01encoder/sampling_3/strided_slice/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2"
 encoder/sampling_3/strided_slice�
encoder/sampling_3/Shape_1Shapeencoder/z_mean/BiasAdd:output:0*
T0*
_output_shapes
:2
encoder/sampling_3/Shape_1�
(encoder/sampling_3/strided_slice_1/stackConst*
_output_shapes
:*
dtype0*
valueB:2*
(encoder/sampling_3/strided_slice_1/stack�
*encoder/sampling_3/strided_slice_1/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2,
*encoder/sampling_3/strided_slice_1/stack_1�
*encoder/sampling_3/strided_slice_1/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2,
*encoder/sampling_3/strided_slice_1/stack_2�
"encoder/sampling_3/strided_slice_1StridedSlice#encoder/sampling_3/Shape_1:output:01encoder/sampling_3/strided_slice_1/stack:output:03encoder/sampling_3/strided_slice_1/stack_1:output:03encoder/sampling_3/strided_slice_1/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2$
"encoder/sampling_3/strided_slice_1�
&encoder/sampling_3/random_normal/shapePack)encoder/sampling_3/strided_slice:output:0+encoder/sampling_3/strided_slice_1:output:0*
N*
T0*
_output_shapes
:2(
&encoder/sampling_3/random_normal/shape�
%encoder/sampling_3/random_normal/meanConst*
_output_shapes
: *
dtype0*
valueB
 *    2'
%encoder/sampling_3/random_normal/mean�
'encoder/sampling_3/random_normal/stddevConst*
_output_shapes
: *
dtype0*
valueB
 *  �?2)
'encoder/sampling_3/random_normal/stddev�
5encoder/sampling_3/random_normal/RandomStandardNormalRandomStandardNormal/encoder/sampling_3/random_normal/shape:output:0*
T0*0
_output_shapes
:������������������*
dtype0*
seed���)*
seed2�;27
5encoder/sampling_3/random_normal/RandomStandardNormal�
$encoder/sampling_3/random_normal/mulMul>encoder/sampling_3/random_normal/RandomStandardNormal:output:00encoder/sampling_3/random_normal/stddev:output:0*
T0*0
_output_shapes
:������������������2&
$encoder/sampling_3/random_normal/mul�
 encoder/sampling_3/random_normalAdd(encoder/sampling_3/random_normal/mul:z:0.encoder/sampling_3/random_normal/mean:output:0*
T0*0
_output_shapes
:������������������2"
 encoder/sampling_3/random_normaly
encoder/sampling_3/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *   ?2
encoder/sampling_3/mul/x�
encoder/sampling_3/mulMul!encoder/sampling_3/mul/x:output:0"encoder/z_log_var/BiasAdd:output:0*
T0*'
_output_shapes
:���������2
encoder/sampling_3/mul�
encoder/sampling_3/ExpExpencoder/sampling_3/mul:z:0*
T0*'
_output_shapes
:���������2
encoder/sampling_3/Exp�
encoder/sampling_3/mul_1Mulencoder/sampling_3/Exp:y:0$encoder/sampling_3/random_normal:z:0*
T0*'
_output_shapes
:���������2
encoder/sampling_3/mul_1�
encoder/sampling_3/addAddV2encoder/z_mean/BiasAdd:output:0encoder/sampling_3/mul_1:z:0*
T0*'
_output_shapes
:���������2
encoder/sampling_3/add{
strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB"    ����2
strided_slice/stack
strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB"        2
strided_slice/stack_1
strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      2
strided_slice/stack_2�
strided_sliceStridedSliceinputsstrided_slice/stack:output:0strided_slice/stack_1:output:0strided_slice/stack_2:output:0*
Index0*
T0*'
_output_shapes
:���������*

begin_mask*
end_mask2
strided_slicet
concatenate/concat/axisConst*
_output_shapes
: *
dtype0*
value	B :2
concatenate/concat/axis�
concatenate/concatConcatV2encoder/sampling_3/add:z:0strided_slice:output:0 concatenate/concat/axis:output:0*
N*
T0*'
_output_shapes
:���������2
concatenate/concat�
(decoder/decoder_l2/MatMul/ReadVariableOpReadVariableOp1decoder_decoder_l2_matmul_readvariableop_resource*
_output_shapes

:*
dtype02*
(decoder/decoder_l2/MatMul/ReadVariableOp�
decoder/decoder_l2/MatMulMatMulconcatenate/concat:output:00decoder/decoder_l2/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2
decoder/decoder_l2/MatMul�
)decoder/decoder_l2/BiasAdd/ReadVariableOpReadVariableOp2decoder_decoder_l2_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02+
)decoder/decoder_l2/BiasAdd/ReadVariableOp�
decoder/decoder_l2/BiasAddBiasAdd#decoder/decoder_l2/MatMul:product:01decoder/decoder_l2/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2
decoder/decoder_l2/BiasAdd�
decoder/decoder_l2/ReluRelu#decoder/decoder_l2/BiasAdd:output:0*
T0*'
_output_shapes
:���������2
decoder/decoder_l2/Relu�
(decoder/decoder_l3/MatMul/ReadVariableOpReadVariableOp1decoder_decoder_l3_matmul_readvariableop_resource*
_output_shapes

:*
dtype02*
(decoder/decoder_l3/MatMul/ReadVariableOp�
decoder/decoder_l3/MatMulMatMul%decoder/decoder_l2/Relu:activations:00decoder/decoder_l3/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2
decoder/decoder_l3/MatMul�
)decoder/decoder_l3/BiasAdd/ReadVariableOpReadVariableOp2decoder_decoder_l3_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02+
)decoder/decoder_l3/BiasAdd/ReadVariableOp�
decoder/decoder_l3/BiasAddBiasAdd#decoder/decoder_l3/MatMul:product:01decoder/decoder_l3/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2
decoder/decoder_l3/BiasAdd�
decoder/decoder_l3/ReluRelu#decoder/decoder_l3/BiasAdd:output:0*
T0*'
_output_shapes
:���������2
decoder/decoder_l3/Relu�
(decoder/decoder_l4/MatMul/ReadVariableOpReadVariableOp1decoder_decoder_l4_matmul_readvariableop_resource*
_output_shapes

:
*
dtype02*
(decoder/decoder_l4/MatMul/ReadVariableOp�
decoder/decoder_l4/MatMulMatMul%decoder/decoder_l3/Relu:activations:00decoder/decoder_l4/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������
2
decoder/decoder_l4/MatMul�
)decoder/decoder_l4/BiasAdd/ReadVariableOpReadVariableOp2decoder_decoder_l4_biasadd_readvariableop_resource*
_output_shapes
:
*
dtype02+
)decoder/decoder_l4/BiasAdd/ReadVariableOp�
decoder/decoder_l4/BiasAddBiasAdd#decoder/decoder_l4/MatMul:product:01decoder/decoder_l4/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������
2
decoder/decoder_l4/BiasAdd�
decoder/decoder_l4/ReluRelu#decoder/decoder_l4/BiasAdd:output:0*
T0*'
_output_shapes
:���������
2
decoder/decoder_l4/Relu�
*decoder/output_layer/MatMul/ReadVariableOpReadVariableOp3decoder_output_layer_matmul_readvariableop_resource*
_output_shapes

:
*
dtype02,
*decoder/output_layer/MatMul/ReadVariableOp�
decoder/output_layer/MatMulMatMul%decoder/decoder_l4/Relu:activations:02decoder/output_layer/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2
decoder/output_layer/MatMul�
+decoder/output_layer/BiasAdd/ReadVariableOpReadVariableOp4decoder_output_layer_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02-
+decoder/output_layer/BiasAdd/ReadVariableOp�
decoder/output_layer/BiasAddBiasAdd%decoder/output_layer/MatMul:product:03decoder/output_layer/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2
decoder/output_layer/BiasAdd�
decoder/output_layer/ReluRelu%decoder/output_layer/BiasAdd:output:0*
T0*'
_output_shapes
:���������2
decoder/output_layer/Relu�
IdentityIdentity'decoder/output_layer/Relu:activations:0*^decoder/decoder_l2/BiasAdd/ReadVariableOp)^decoder/decoder_l2/MatMul/ReadVariableOp*^decoder/decoder_l3/BiasAdd/ReadVariableOp)^decoder/decoder_l3/MatMul/ReadVariableOp*^decoder/decoder_l4/BiasAdd/ReadVariableOp)^decoder/decoder_l4/MatMul/ReadVariableOp,^decoder/output_layer/BiasAdd/ReadVariableOp+^decoder/output_layer/MatMul/ReadVariableOp*^encoder/encoder_l2/BiasAdd/ReadVariableOp)^encoder/encoder_l2/MatMul/ReadVariableOp*^encoder/encoder_l3/BiasAdd/ReadVariableOp)^encoder/encoder_l3/MatMul/ReadVariableOp*^encoder/encoder_l4/BiasAdd/ReadVariableOp)^encoder/encoder_l4/MatMul/ReadVariableOp)^encoder/z_log_var/BiasAdd/ReadVariableOp(^encoder/z_log_var/MatMul/ReadVariableOp&^encoder/z_mean/BiasAdd/ReadVariableOp%^encoder/z_mean/MatMul/ReadVariableOp*
T0*'
_output_shapes
:���������2

Identity"
identityIdentity:output:0*n
_input_shapes]
[:���������::::::::::::::::::2V
)decoder/decoder_l2/BiasAdd/ReadVariableOp)decoder/decoder_l2/BiasAdd/ReadVariableOp2T
(decoder/decoder_l2/MatMul/ReadVariableOp(decoder/decoder_l2/MatMul/ReadVariableOp2V
)decoder/decoder_l3/BiasAdd/ReadVariableOp)decoder/decoder_l3/BiasAdd/ReadVariableOp2T
(decoder/decoder_l3/MatMul/ReadVariableOp(decoder/decoder_l3/MatMul/ReadVariableOp2V
)decoder/decoder_l4/BiasAdd/ReadVariableOp)decoder/decoder_l4/BiasAdd/ReadVariableOp2T
(decoder/decoder_l4/MatMul/ReadVariableOp(decoder/decoder_l4/MatMul/ReadVariableOp2Z
+decoder/output_layer/BiasAdd/ReadVariableOp+decoder/output_layer/BiasAdd/ReadVariableOp2X
*decoder/output_layer/MatMul/ReadVariableOp*decoder/output_layer/MatMul/ReadVariableOp2V
)encoder/encoder_l2/BiasAdd/ReadVariableOp)encoder/encoder_l2/BiasAdd/ReadVariableOp2T
(encoder/encoder_l2/MatMul/ReadVariableOp(encoder/encoder_l2/MatMul/ReadVariableOp2V
)encoder/encoder_l3/BiasAdd/ReadVariableOp)encoder/encoder_l3/BiasAdd/ReadVariableOp2T
(encoder/encoder_l3/MatMul/ReadVariableOp(encoder/encoder_l3/MatMul/ReadVariableOp2V
)encoder/encoder_l4/BiasAdd/ReadVariableOp)encoder/encoder_l4/BiasAdd/ReadVariableOp2T
(encoder/encoder_l4/MatMul/ReadVariableOp(encoder/encoder_l4/MatMul/ReadVariableOp2T
(encoder/z_log_var/BiasAdd/ReadVariableOp(encoder/z_log_var/BiasAdd/ReadVariableOp2R
'encoder/z_log_var/MatMul/ReadVariableOp'encoder/z_log_var/MatMul/ReadVariableOp2N
%encoder/z_mean/BiasAdd/ReadVariableOp%encoder/z_mean/BiasAdd/ReadVariableOp2L
$encoder/z_mean/MatMul/ReadVariableOp$encoder/z_mean/MatMul/ReadVariableOp:O K
'
_output_shapes
:���������
 
_user_specified_nameinputs
�	
�
F__inference_encoder_l3_layer_call_and_return_conditional_losses_141194

inputs"
matmul_readvariableop_resource#
biasadd_readvariableop_resource
identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOp�
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:
*
dtype02
MatMul/ReadVariableOps
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2
MatMul�
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:*
dtype02
BiasAdd/ReadVariableOp�
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2	
BiasAddX
ReluReluBiasAdd:output:0*
T0*'
_output_shapes
:���������2
Relu�
IdentityIdentityRelu:activations:0^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*
T0*'
_output_shapes
:���������2

Identity"
identityIdentity:output:0*.
_input_shapes
:���������
::20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:O K
'
_output_shapes
:���������

 
_user_specified_nameinputs
�	
�
B__inference_z_mean_layer_call_and_return_conditional_losses_142981

inputs"
matmul_readvariableop_resource#
biasadd_readvariableop_resource
identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOp�
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:*
dtype02
MatMul/ReadVariableOps
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2
MatMul�
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:*
dtype02
BiasAdd/ReadVariableOp�
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2	
BiasAdd�
IdentityIdentityBiasAdd:output:0^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*
T0*'
_output_shapes
:���������2

Identity"
identityIdentity:output:0*.
_input_shapes
:���������::20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:O K
'
_output_shapes
:���������
 
_user_specified_nameinputs
�&
�
C__inference_encoder_layer_call_and_return_conditional_losses_141359
encoder_input
encoder_l2_141330
encoder_l2_141332
encoder_l3_141335
encoder_l3_141337
encoder_l4_141340
encoder_l4_141342
z_mean_141345
z_mean_141347
z_log_var_141350
z_log_var_141352
identity

identity_1

identity_2��"encoder_l2/StatefulPartitionedCall�"encoder_l3/StatefulPartitionedCall�"encoder_l4/StatefulPartitionedCall�"sampling_3/StatefulPartitionedCall�!z_log_var/StatefulPartitionedCall�z_mean/StatefulPartitionedCall�
"encoder_l2/StatefulPartitionedCallStatefulPartitionedCallencoder_inputencoder_l2_141330encoder_l2_141332*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������
*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *O
fJRH
F__inference_encoder_l2_layer_call_and_return_conditional_losses_1411672$
"encoder_l2/StatefulPartitionedCall�
"encoder_l3/StatefulPartitionedCallStatefulPartitionedCall+encoder_l2/StatefulPartitionedCall:output:0encoder_l3_141335encoder_l3_141337*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *O
fJRH
F__inference_encoder_l3_layer_call_and_return_conditional_losses_1411942$
"encoder_l3/StatefulPartitionedCall�
"encoder_l4/StatefulPartitionedCallStatefulPartitionedCall+encoder_l3/StatefulPartitionedCall:output:0encoder_l4_141340encoder_l4_141342*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *O
fJRH
F__inference_encoder_l4_layer_call_and_return_conditional_losses_1412212$
"encoder_l4/StatefulPartitionedCall�
z_mean/StatefulPartitionedCallStatefulPartitionedCall+encoder_l4/StatefulPartitionedCall:output:0z_mean_141345z_mean_141347*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *K
fFRD
B__inference_z_mean_layer_call_and_return_conditional_losses_1412472 
z_mean/StatefulPartitionedCall�
!z_log_var/StatefulPartitionedCallStatefulPartitionedCall+encoder_l4/StatefulPartitionedCall:output:0z_log_var_141350z_log_var_141352*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *N
fIRG
E__inference_z_log_var_layer_call_and_return_conditional_losses_1412732#
!z_log_var/StatefulPartitionedCall�
"sampling_3/StatefulPartitionedCallStatefulPartitionedCall'z_mean/StatefulPartitionedCall:output:0*z_log_var/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8� *O
fJRH
F__inference_sampling_3_layer_call_and_return_conditional_losses_1413152$
"sampling_3/StatefulPartitionedCall�
IdentityIdentity'z_mean/StatefulPartitionedCall:output:0#^encoder_l2/StatefulPartitionedCall#^encoder_l3/StatefulPartitionedCall#^encoder_l4/StatefulPartitionedCall#^sampling_3/StatefulPartitionedCall"^z_log_var/StatefulPartitionedCall^z_mean/StatefulPartitionedCall*
T0*'
_output_shapes
:���������2

Identity�

Identity_1Identity*z_log_var/StatefulPartitionedCall:output:0#^encoder_l2/StatefulPartitionedCall#^encoder_l3/StatefulPartitionedCall#^encoder_l4/StatefulPartitionedCall#^sampling_3/StatefulPartitionedCall"^z_log_var/StatefulPartitionedCall^z_mean/StatefulPartitionedCall*
T0*'
_output_shapes
:���������2

Identity_1�

Identity_2Identity+sampling_3/StatefulPartitionedCall:output:0#^encoder_l2/StatefulPartitionedCall#^encoder_l3/StatefulPartitionedCall#^encoder_l4/StatefulPartitionedCall#^sampling_3/StatefulPartitionedCall"^z_log_var/StatefulPartitionedCall^z_mean/StatefulPartitionedCall*
T0*'
_output_shapes
:���������2

Identity_2"
identityIdentity:output:0"!

identity_1Identity_1:output:0"!

identity_2Identity_2:output:0*N
_input_shapes=
;:���������::::::::::2H
"encoder_l2/StatefulPartitionedCall"encoder_l2/StatefulPartitionedCall2H
"encoder_l3/StatefulPartitionedCall"encoder_l3/StatefulPartitionedCall2H
"encoder_l4/StatefulPartitionedCall"encoder_l4/StatefulPartitionedCall2H
"sampling_3/StatefulPartitionedCall"sampling_3/StatefulPartitionedCall2F
!z_log_var/StatefulPartitionedCall!z_log_var/StatefulPartitionedCall2@
z_mean/StatefulPartitionedCallz_mean/StatefulPartitionedCall:V R
'
_output_shapes
:���������
'
_user_specified_nameencoder_input
�
�
(__inference_decoder_layer_call_fn_142911

inputs
unknown
	unknown_0
	unknown_1
	unknown_2
	unknown_3
	unknown_4
	unknown_5
	unknown_6
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6*
Tin
2	*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������**
_read_only_resource_inputs

*-
config_proto

CPU

GPU 2J 8� *L
fGRE
C__inference_decoder_layer_call_and_return_conditional_losses_1416912
StatefulPartitionedCall�
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*'
_output_shapes
:���������2

Identity"
identityIdentity:output:0*F
_input_shapes5
3:���������::::::::22
StatefulPartitionedCallStatefulPartitionedCall:O K
'
_output_shapes
:���������
 
_user_specified_nameinputs
�	
�
F__inference_decoder_l3_layer_call_and_return_conditional_losses_141524

inputs"
matmul_readvariableop_resource#
biasadd_readvariableop_resource
identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOp�
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:*
dtype02
MatMul/ReadVariableOps
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2
MatMul�
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:*
dtype02
BiasAdd/ReadVariableOp�
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2	
BiasAddX
ReluReluBiasAdd:output:0*
T0*'
_output_shapes
:���������2
Relu�
IdentityIdentityRelu:activations:0^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*
T0*'
_output_shapes
:���������2

Identity"
identityIdentity:output:0*.
_input_shapes
:���������::20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:O K
'
_output_shapes
:���������
 
_user_specified_nameinputs
�
�
C__inference_decoder_layer_call_and_return_conditional_losses_141646

inputs
decoder_l2_141625
decoder_l2_141627
decoder_l3_141630
decoder_l3_141632
decoder_l4_141635
decoder_l4_141637
output_layer_141640
output_layer_141642
identity��"decoder_l2/StatefulPartitionedCall�"decoder_l3/StatefulPartitionedCall�"decoder_l4/StatefulPartitionedCall�$output_layer/StatefulPartitionedCall�
"decoder_l2/StatefulPartitionedCallStatefulPartitionedCallinputsdecoder_l2_141625decoder_l2_141627*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *O
fJRH
F__inference_decoder_l2_layer_call_and_return_conditional_losses_1414972$
"decoder_l2/StatefulPartitionedCall�
"decoder_l3/StatefulPartitionedCallStatefulPartitionedCall+decoder_l2/StatefulPartitionedCall:output:0decoder_l3_141630decoder_l3_141632*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *O
fJRH
F__inference_decoder_l3_layer_call_and_return_conditional_losses_1415242$
"decoder_l3/StatefulPartitionedCall�
"decoder_l4/StatefulPartitionedCallStatefulPartitionedCall+decoder_l3/StatefulPartitionedCall:output:0decoder_l4_141635decoder_l4_141637*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������
*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *O
fJRH
F__inference_decoder_l4_layer_call_and_return_conditional_losses_1415512$
"decoder_l4/StatefulPartitionedCall�
$output_layer/StatefulPartitionedCallStatefulPartitionedCall+decoder_l4/StatefulPartitionedCall:output:0output_layer_141640output_layer_141642*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *Q
fLRJ
H__inference_output_layer_layer_call_and_return_conditional_losses_1415782&
$output_layer/StatefulPartitionedCall�
IdentityIdentity-output_layer/StatefulPartitionedCall:output:0#^decoder_l2/StatefulPartitionedCall#^decoder_l3/StatefulPartitionedCall#^decoder_l4/StatefulPartitionedCall%^output_layer/StatefulPartitionedCall*
T0*'
_output_shapes
:���������2

Identity"
identityIdentity:output:0*F
_input_shapes5
3:���������::::::::2H
"decoder_l2/StatefulPartitionedCall"decoder_l2/StatefulPartitionedCall2H
"decoder_l3/StatefulPartitionedCall"decoder_l3/StatefulPartitionedCall2H
"decoder_l4/StatefulPartitionedCall"decoder_l4/StatefulPartitionedCall2L
$output_layer/StatefulPartitionedCall$output_layer/StatefulPartitionedCall:O K
'
_output_shapes
:���������
 
_user_specified_nameinputs
��
�
B__inference_cvae_3_layer_call_and_return_conditional_losses_142187
input_15
1encoder_encoder_l2_matmul_readvariableop_resource6
2encoder_encoder_l2_biasadd_readvariableop_resource5
1encoder_encoder_l3_matmul_readvariableop_resource6
2encoder_encoder_l3_biasadd_readvariableop_resource5
1encoder_encoder_l4_matmul_readvariableop_resource6
2encoder_encoder_l4_biasadd_readvariableop_resource1
-encoder_z_mean_matmul_readvariableop_resource2
.encoder_z_mean_biasadd_readvariableop_resource4
0encoder_z_log_var_matmul_readvariableop_resource5
1encoder_z_log_var_biasadd_readvariableop_resource5
1decoder_decoder_l2_matmul_readvariableop_resource6
2decoder_decoder_l2_biasadd_readvariableop_resource5
1decoder_decoder_l3_matmul_readvariableop_resource6
2decoder_decoder_l3_biasadd_readvariableop_resource5
1decoder_decoder_l4_matmul_readvariableop_resource6
2decoder_decoder_l4_biasadd_readvariableop_resource7
3decoder_output_layer_matmul_readvariableop_resource8
4decoder_output_layer_biasadd_readvariableop_resource
identity��)decoder/decoder_l2/BiasAdd/ReadVariableOp�(decoder/decoder_l2/MatMul/ReadVariableOp�)decoder/decoder_l3/BiasAdd/ReadVariableOp�(decoder/decoder_l3/MatMul/ReadVariableOp�)decoder/decoder_l4/BiasAdd/ReadVariableOp�(decoder/decoder_l4/MatMul/ReadVariableOp�+decoder/output_layer/BiasAdd/ReadVariableOp�*decoder/output_layer/MatMul/ReadVariableOp�)encoder/encoder_l2/BiasAdd/ReadVariableOp�(encoder/encoder_l2/MatMul/ReadVariableOp�)encoder/encoder_l3/BiasAdd/ReadVariableOp�(encoder/encoder_l3/MatMul/ReadVariableOp�)encoder/encoder_l4/BiasAdd/ReadVariableOp�(encoder/encoder_l4/MatMul/ReadVariableOp�(encoder/z_log_var/BiasAdd/ReadVariableOp�'encoder/z_log_var/MatMul/ReadVariableOp�%encoder/z_mean/BiasAdd/ReadVariableOp�$encoder/z_mean/MatMul/ReadVariableOp�
(encoder/encoder_l2/MatMul/ReadVariableOpReadVariableOp1encoder_encoder_l2_matmul_readvariableop_resource*
_output_shapes

:
*
dtype02*
(encoder/encoder_l2/MatMul/ReadVariableOp�
encoder/encoder_l2/MatMulMatMulinput_10encoder/encoder_l2/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������
2
encoder/encoder_l2/MatMul�
)encoder/encoder_l2/BiasAdd/ReadVariableOpReadVariableOp2encoder_encoder_l2_biasadd_readvariableop_resource*
_output_shapes
:
*
dtype02+
)encoder/encoder_l2/BiasAdd/ReadVariableOp�
encoder/encoder_l2/BiasAddBiasAdd#encoder/encoder_l2/MatMul:product:01encoder/encoder_l2/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������
2
encoder/encoder_l2/BiasAdd�
encoder/encoder_l2/ReluRelu#encoder/encoder_l2/BiasAdd:output:0*
T0*'
_output_shapes
:���������
2
encoder/encoder_l2/Relu�
(encoder/encoder_l3/MatMul/ReadVariableOpReadVariableOp1encoder_encoder_l3_matmul_readvariableop_resource*
_output_shapes

:
*
dtype02*
(encoder/encoder_l3/MatMul/ReadVariableOp�
encoder/encoder_l3/MatMulMatMul%encoder/encoder_l2/Relu:activations:00encoder/encoder_l3/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2
encoder/encoder_l3/MatMul�
)encoder/encoder_l3/BiasAdd/ReadVariableOpReadVariableOp2encoder_encoder_l3_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02+
)encoder/encoder_l3/BiasAdd/ReadVariableOp�
encoder/encoder_l3/BiasAddBiasAdd#encoder/encoder_l3/MatMul:product:01encoder/encoder_l3/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2
encoder/encoder_l3/BiasAdd�
encoder/encoder_l3/ReluRelu#encoder/encoder_l3/BiasAdd:output:0*
T0*'
_output_shapes
:���������2
encoder/encoder_l3/Relu�
(encoder/encoder_l4/MatMul/ReadVariableOpReadVariableOp1encoder_encoder_l4_matmul_readvariableop_resource*
_output_shapes

:*
dtype02*
(encoder/encoder_l4/MatMul/ReadVariableOp�
encoder/encoder_l4/MatMulMatMul%encoder/encoder_l3/Relu:activations:00encoder/encoder_l4/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2
encoder/encoder_l4/MatMul�
)encoder/encoder_l4/BiasAdd/ReadVariableOpReadVariableOp2encoder_encoder_l4_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02+
)encoder/encoder_l4/BiasAdd/ReadVariableOp�
encoder/encoder_l4/BiasAddBiasAdd#encoder/encoder_l4/MatMul:product:01encoder/encoder_l4/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2
encoder/encoder_l4/BiasAdd�
encoder/encoder_l4/ReluRelu#encoder/encoder_l4/BiasAdd:output:0*
T0*'
_output_shapes
:���������2
encoder/encoder_l4/Relu�
$encoder/z_mean/MatMul/ReadVariableOpReadVariableOp-encoder_z_mean_matmul_readvariableop_resource*
_output_shapes

:*
dtype02&
$encoder/z_mean/MatMul/ReadVariableOp�
encoder/z_mean/MatMulMatMul%encoder/encoder_l4/Relu:activations:0,encoder/z_mean/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2
encoder/z_mean/MatMul�
%encoder/z_mean/BiasAdd/ReadVariableOpReadVariableOp.encoder_z_mean_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02'
%encoder/z_mean/BiasAdd/ReadVariableOp�
encoder/z_mean/BiasAddBiasAddencoder/z_mean/MatMul:product:0-encoder/z_mean/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2
encoder/z_mean/BiasAdd�
'encoder/z_log_var/MatMul/ReadVariableOpReadVariableOp0encoder_z_log_var_matmul_readvariableop_resource*
_output_shapes

:*
dtype02)
'encoder/z_log_var/MatMul/ReadVariableOp�
encoder/z_log_var/MatMulMatMul%encoder/encoder_l4/Relu:activations:0/encoder/z_log_var/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2
encoder/z_log_var/MatMul�
(encoder/z_log_var/BiasAdd/ReadVariableOpReadVariableOp1encoder_z_log_var_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02*
(encoder/z_log_var/BiasAdd/ReadVariableOp�
encoder/z_log_var/BiasAddBiasAdd"encoder/z_log_var/MatMul:product:00encoder/z_log_var/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2
encoder/z_log_var/BiasAdd�
encoder/sampling_3/ShapeShapeencoder/z_mean/BiasAdd:output:0*
T0*
_output_shapes
:2
encoder/sampling_3/Shape�
&encoder/sampling_3/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB: 2(
&encoder/sampling_3/strided_slice/stack�
(encoder/sampling_3/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2*
(encoder/sampling_3/strided_slice/stack_1�
(encoder/sampling_3/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2*
(encoder/sampling_3/strided_slice/stack_2�
 encoder/sampling_3/strided_sliceStridedSlice!encoder/sampling_3/Shape:output:0/encoder/sampling_3/strided_slice/stack:output:01encoder/sampling_3/strided_slice/stack_1:output:01encoder/sampling_3/strided_slice/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2"
 encoder/sampling_3/strided_slice�
encoder/sampling_3/Shape_1Shapeencoder/z_mean/BiasAdd:output:0*
T0*
_output_shapes
:2
encoder/sampling_3/Shape_1�
(encoder/sampling_3/strided_slice_1/stackConst*
_output_shapes
:*
dtype0*
valueB:2*
(encoder/sampling_3/strided_slice_1/stack�
*encoder/sampling_3/strided_slice_1/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2,
*encoder/sampling_3/strided_slice_1/stack_1�
*encoder/sampling_3/strided_slice_1/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2,
*encoder/sampling_3/strided_slice_1/stack_2�
"encoder/sampling_3/strided_slice_1StridedSlice#encoder/sampling_3/Shape_1:output:01encoder/sampling_3/strided_slice_1/stack:output:03encoder/sampling_3/strided_slice_1/stack_1:output:03encoder/sampling_3/strided_slice_1/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2$
"encoder/sampling_3/strided_slice_1�
&encoder/sampling_3/random_normal/shapePack)encoder/sampling_3/strided_slice:output:0+encoder/sampling_3/strided_slice_1:output:0*
N*
T0*
_output_shapes
:2(
&encoder/sampling_3/random_normal/shape�
%encoder/sampling_3/random_normal/meanConst*
_output_shapes
: *
dtype0*
valueB
 *    2'
%encoder/sampling_3/random_normal/mean�
'encoder/sampling_3/random_normal/stddevConst*
_output_shapes
: *
dtype0*
valueB
 *  �?2)
'encoder/sampling_3/random_normal/stddev�
5encoder/sampling_3/random_normal/RandomStandardNormalRandomStandardNormal/encoder/sampling_3/random_normal/shape:output:0*
T0*0
_output_shapes
:������������������*
dtype0*
seed���)*
seed2��27
5encoder/sampling_3/random_normal/RandomStandardNormal�
$encoder/sampling_3/random_normal/mulMul>encoder/sampling_3/random_normal/RandomStandardNormal:output:00encoder/sampling_3/random_normal/stddev:output:0*
T0*0
_output_shapes
:������������������2&
$encoder/sampling_3/random_normal/mul�
 encoder/sampling_3/random_normalAdd(encoder/sampling_3/random_normal/mul:z:0.encoder/sampling_3/random_normal/mean:output:0*
T0*0
_output_shapes
:������������������2"
 encoder/sampling_3/random_normaly
encoder/sampling_3/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *   ?2
encoder/sampling_3/mul/x�
encoder/sampling_3/mulMul!encoder/sampling_3/mul/x:output:0"encoder/z_log_var/BiasAdd:output:0*
T0*'
_output_shapes
:���������2
encoder/sampling_3/mul�
encoder/sampling_3/ExpExpencoder/sampling_3/mul:z:0*
T0*'
_output_shapes
:���������2
encoder/sampling_3/Exp�
encoder/sampling_3/mul_1Mulencoder/sampling_3/Exp:y:0$encoder/sampling_3/random_normal:z:0*
T0*'
_output_shapes
:���������2
encoder/sampling_3/mul_1�
encoder/sampling_3/addAddV2encoder/z_mean/BiasAdd:output:0encoder/sampling_3/mul_1:z:0*
T0*'
_output_shapes
:���������2
encoder/sampling_3/add{
strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB"    ����2
strided_slice/stack
strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB"        2
strided_slice/stack_1
strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      2
strided_slice/stack_2�
strided_sliceStridedSliceinput_1strided_slice/stack:output:0strided_slice/stack_1:output:0strided_slice/stack_2:output:0*
Index0*
T0*'
_output_shapes
:���������*

begin_mask*
end_mask2
strided_slicet
concatenate/concat/axisConst*
_output_shapes
: *
dtype0*
value	B :2
concatenate/concat/axis�
concatenate/concatConcatV2encoder/sampling_3/add:z:0strided_slice:output:0 concatenate/concat/axis:output:0*
N*
T0*'
_output_shapes
:���������2
concatenate/concat�
(decoder/decoder_l2/MatMul/ReadVariableOpReadVariableOp1decoder_decoder_l2_matmul_readvariableop_resource*
_output_shapes

:*
dtype02*
(decoder/decoder_l2/MatMul/ReadVariableOp�
decoder/decoder_l2/MatMulMatMulconcatenate/concat:output:00decoder/decoder_l2/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2
decoder/decoder_l2/MatMul�
)decoder/decoder_l2/BiasAdd/ReadVariableOpReadVariableOp2decoder_decoder_l2_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02+
)decoder/decoder_l2/BiasAdd/ReadVariableOp�
decoder/decoder_l2/BiasAddBiasAdd#decoder/decoder_l2/MatMul:product:01decoder/decoder_l2/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2
decoder/decoder_l2/BiasAdd�
decoder/decoder_l2/ReluRelu#decoder/decoder_l2/BiasAdd:output:0*
T0*'
_output_shapes
:���������2
decoder/decoder_l2/Relu�
(decoder/decoder_l3/MatMul/ReadVariableOpReadVariableOp1decoder_decoder_l3_matmul_readvariableop_resource*
_output_shapes

:*
dtype02*
(decoder/decoder_l3/MatMul/ReadVariableOp�
decoder/decoder_l3/MatMulMatMul%decoder/decoder_l2/Relu:activations:00decoder/decoder_l3/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2
decoder/decoder_l3/MatMul�
)decoder/decoder_l3/BiasAdd/ReadVariableOpReadVariableOp2decoder_decoder_l3_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02+
)decoder/decoder_l3/BiasAdd/ReadVariableOp�
decoder/decoder_l3/BiasAddBiasAdd#decoder/decoder_l3/MatMul:product:01decoder/decoder_l3/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2
decoder/decoder_l3/BiasAdd�
decoder/decoder_l3/ReluRelu#decoder/decoder_l3/BiasAdd:output:0*
T0*'
_output_shapes
:���������2
decoder/decoder_l3/Relu�
(decoder/decoder_l4/MatMul/ReadVariableOpReadVariableOp1decoder_decoder_l4_matmul_readvariableop_resource*
_output_shapes

:
*
dtype02*
(decoder/decoder_l4/MatMul/ReadVariableOp�
decoder/decoder_l4/MatMulMatMul%decoder/decoder_l3/Relu:activations:00decoder/decoder_l4/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������
2
decoder/decoder_l4/MatMul�
)decoder/decoder_l4/BiasAdd/ReadVariableOpReadVariableOp2decoder_decoder_l4_biasadd_readvariableop_resource*
_output_shapes
:
*
dtype02+
)decoder/decoder_l4/BiasAdd/ReadVariableOp�
decoder/decoder_l4/BiasAddBiasAdd#decoder/decoder_l4/MatMul:product:01decoder/decoder_l4/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������
2
decoder/decoder_l4/BiasAdd�
decoder/decoder_l4/ReluRelu#decoder/decoder_l4/BiasAdd:output:0*
T0*'
_output_shapes
:���������
2
decoder/decoder_l4/Relu�
*decoder/output_layer/MatMul/ReadVariableOpReadVariableOp3decoder_output_layer_matmul_readvariableop_resource*
_output_shapes

:
*
dtype02,
*decoder/output_layer/MatMul/ReadVariableOp�
decoder/output_layer/MatMulMatMul%decoder/decoder_l4/Relu:activations:02decoder/output_layer/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2
decoder/output_layer/MatMul�
+decoder/output_layer/BiasAdd/ReadVariableOpReadVariableOp4decoder_output_layer_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02-
+decoder/output_layer/BiasAdd/ReadVariableOp�
decoder/output_layer/BiasAddBiasAdd%decoder/output_layer/MatMul:product:03decoder/output_layer/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2
decoder/output_layer/BiasAdd�
decoder/output_layer/ReluRelu%decoder/output_layer/BiasAdd:output:0*
T0*'
_output_shapes
:���������2
decoder/output_layer/Relu�
IdentityIdentity'decoder/output_layer/Relu:activations:0*^decoder/decoder_l2/BiasAdd/ReadVariableOp)^decoder/decoder_l2/MatMul/ReadVariableOp*^decoder/decoder_l3/BiasAdd/ReadVariableOp)^decoder/decoder_l3/MatMul/ReadVariableOp*^decoder/decoder_l4/BiasAdd/ReadVariableOp)^decoder/decoder_l4/MatMul/ReadVariableOp,^decoder/output_layer/BiasAdd/ReadVariableOp+^decoder/output_layer/MatMul/ReadVariableOp*^encoder/encoder_l2/BiasAdd/ReadVariableOp)^encoder/encoder_l2/MatMul/ReadVariableOp*^encoder/encoder_l3/BiasAdd/ReadVariableOp)^encoder/encoder_l3/MatMul/ReadVariableOp*^encoder/encoder_l4/BiasAdd/ReadVariableOp)^encoder/encoder_l4/MatMul/ReadVariableOp)^encoder/z_log_var/BiasAdd/ReadVariableOp(^encoder/z_log_var/MatMul/ReadVariableOp&^encoder/z_mean/BiasAdd/ReadVariableOp%^encoder/z_mean/MatMul/ReadVariableOp*
T0*'
_output_shapes
:���������2

Identity"
identityIdentity:output:0*n
_input_shapes]
[:���������::::::::::::::::::2V
)decoder/decoder_l2/BiasAdd/ReadVariableOp)decoder/decoder_l2/BiasAdd/ReadVariableOp2T
(decoder/decoder_l2/MatMul/ReadVariableOp(decoder/decoder_l2/MatMul/ReadVariableOp2V
)decoder/decoder_l3/BiasAdd/ReadVariableOp)decoder/decoder_l3/BiasAdd/ReadVariableOp2T
(decoder/decoder_l3/MatMul/ReadVariableOp(decoder/decoder_l3/MatMul/ReadVariableOp2V
)decoder/decoder_l4/BiasAdd/ReadVariableOp)decoder/decoder_l4/BiasAdd/ReadVariableOp2T
(decoder/decoder_l4/MatMul/ReadVariableOp(decoder/decoder_l4/MatMul/ReadVariableOp2Z
+decoder/output_layer/BiasAdd/ReadVariableOp+decoder/output_layer/BiasAdd/ReadVariableOp2X
*decoder/output_layer/MatMul/ReadVariableOp*decoder/output_layer/MatMul/ReadVariableOp2V
)encoder/encoder_l2/BiasAdd/ReadVariableOp)encoder/encoder_l2/BiasAdd/ReadVariableOp2T
(encoder/encoder_l2/MatMul/ReadVariableOp(encoder/encoder_l2/MatMul/ReadVariableOp2V
)encoder/encoder_l3/BiasAdd/ReadVariableOp)encoder/encoder_l3/BiasAdd/ReadVariableOp2T
(encoder/encoder_l3/MatMul/ReadVariableOp(encoder/encoder_l3/MatMul/ReadVariableOp2V
)encoder/encoder_l4/BiasAdd/ReadVariableOp)encoder/encoder_l4/BiasAdd/ReadVariableOp2T
(encoder/encoder_l4/MatMul/ReadVariableOp(encoder/encoder_l4/MatMul/ReadVariableOp2T
(encoder/z_log_var/BiasAdd/ReadVariableOp(encoder/z_log_var/BiasAdd/ReadVariableOp2R
'encoder/z_log_var/MatMul/ReadVariableOp'encoder/z_log_var/MatMul/ReadVariableOp2N
%encoder/z_mean/BiasAdd/ReadVariableOp%encoder/z_mean/BiasAdd/ReadVariableOp2L
$encoder/z_mean/MatMul/ReadVariableOp$encoder/z_mean/MatMul/ReadVariableOp:P L
'
_output_shapes
:���������
!
_user_specified_name	input_1
�
u
F__inference_sampling_3_layer_call_and_return_conditional_losses_143035
inputs_0
inputs_1
identity�F
ShapeShapeinputs_0*
T0*
_output_shapes
:2
Shapet
strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB: 2
strided_slice/stackx
strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice/stack_1x
strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice/stack_2�
strided_sliceStridedSliceShape:output:0strided_slice/stack:output:0strided_slice/stack_1:output:0strided_slice/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2
strided_sliceJ
Shape_1Shapeinputs_0*
T0*
_output_shapes
:2	
Shape_1x
strided_slice_1/stackConst*
_output_shapes
:*
dtype0*
valueB:2
strided_slice_1/stack|
strided_slice_1/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice_1/stack_1|
strided_slice_1/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice_1/stack_2�
strided_slice_1StridedSliceShape_1:output:0strided_slice_1/stack:output:0 strided_slice_1/stack_1:output:0 strided_slice_1/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2
strided_slice_1�
random_normal/shapePackstrided_slice:output:0strided_slice_1:output:0*
N*
T0*
_output_shapes
:2
random_normal/shapem
random_normal/meanConst*
_output_shapes
: *
dtype0*
valueB
 *    2
random_normal/meanq
random_normal/stddevConst*
_output_shapes
: *
dtype0*
valueB
 *  �?2
random_normal/stddev�
"random_normal/RandomStandardNormalRandomStandardNormalrandom_normal/shape:output:0*
T0*0
_output_shapes
:������������������*
dtype0*
seed���)*
seed2���2$
"random_normal/RandomStandardNormal�
random_normal/mulMul+random_normal/RandomStandardNormal:output:0random_normal/stddev:output:0*
T0*0
_output_shapes
:������������������2
random_normal/mul�
random_normalAddrandom_normal/mul:z:0random_normal/mean:output:0*
T0*0
_output_shapes
:������������������2
random_normalS
mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *   ?2
mul/x]
mulMulmul/x:output:0inputs_1*
T0*'
_output_shapes
:���������2
mulL
ExpExpmul:z:0*
T0*'
_output_shapes
:���������2
Expc
mul_1MulExp:y:0random_normal:z:0*
T0*'
_output_shapes
:���������2
mul_1Z
addAddV2inputs_0	mul_1:z:0*
T0*'
_output_shapes
:���������2
add[
IdentityIdentityadd:z:0*
T0*'
_output_shapes
:���������2

Identity"
identityIdentity:output:0*9
_input_shapes(
&:���������:���������:Q M
'
_output_shapes
:���������
"
_user_specified_name
inputs/0:QM
'
_output_shapes
:���������
"
_user_specified_name
inputs/1
՜
�
!__inference__wrapped_model_141152
input_1<
8cvae_3_encoder_encoder_l2_matmul_readvariableop_resource=
9cvae_3_encoder_encoder_l2_biasadd_readvariableop_resource<
8cvae_3_encoder_encoder_l3_matmul_readvariableop_resource=
9cvae_3_encoder_encoder_l3_biasadd_readvariableop_resource<
8cvae_3_encoder_encoder_l4_matmul_readvariableop_resource=
9cvae_3_encoder_encoder_l4_biasadd_readvariableop_resource8
4cvae_3_encoder_z_mean_matmul_readvariableop_resource9
5cvae_3_encoder_z_mean_biasadd_readvariableop_resource;
7cvae_3_encoder_z_log_var_matmul_readvariableop_resource<
8cvae_3_encoder_z_log_var_biasadd_readvariableop_resource<
8cvae_3_decoder_decoder_l2_matmul_readvariableop_resource=
9cvae_3_decoder_decoder_l2_biasadd_readvariableop_resource<
8cvae_3_decoder_decoder_l3_matmul_readvariableop_resource=
9cvae_3_decoder_decoder_l3_biasadd_readvariableop_resource<
8cvae_3_decoder_decoder_l4_matmul_readvariableop_resource=
9cvae_3_decoder_decoder_l4_biasadd_readvariableop_resource>
:cvae_3_decoder_output_layer_matmul_readvariableop_resource?
;cvae_3_decoder_output_layer_biasadd_readvariableop_resource
identity��0cvae_3/decoder/decoder_l2/BiasAdd/ReadVariableOp�/cvae_3/decoder/decoder_l2/MatMul/ReadVariableOp�0cvae_3/decoder/decoder_l3/BiasAdd/ReadVariableOp�/cvae_3/decoder/decoder_l3/MatMul/ReadVariableOp�0cvae_3/decoder/decoder_l4/BiasAdd/ReadVariableOp�/cvae_3/decoder/decoder_l4/MatMul/ReadVariableOp�2cvae_3/decoder/output_layer/BiasAdd/ReadVariableOp�1cvae_3/decoder/output_layer/MatMul/ReadVariableOp�0cvae_3/encoder/encoder_l2/BiasAdd/ReadVariableOp�/cvae_3/encoder/encoder_l2/MatMul/ReadVariableOp�0cvae_3/encoder/encoder_l3/BiasAdd/ReadVariableOp�/cvae_3/encoder/encoder_l3/MatMul/ReadVariableOp�0cvae_3/encoder/encoder_l4/BiasAdd/ReadVariableOp�/cvae_3/encoder/encoder_l4/MatMul/ReadVariableOp�/cvae_3/encoder/z_log_var/BiasAdd/ReadVariableOp�.cvae_3/encoder/z_log_var/MatMul/ReadVariableOp�,cvae_3/encoder/z_mean/BiasAdd/ReadVariableOp�+cvae_3/encoder/z_mean/MatMul/ReadVariableOp�
/cvae_3/encoder/encoder_l2/MatMul/ReadVariableOpReadVariableOp8cvae_3_encoder_encoder_l2_matmul_readvariableop_resource*
_output_shapes

:
*
dtype021
/cvae_3/encoder/encoder_l2/MatMul/ReadVariableOp�
 cvae_3/encoder/encoder_l2/MatMulMatMulinput_17cvae_3/encoder/encoder_l2/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������
2"
 cvae_3/encoder/encoder_l2/MatMul�
0cvae_3/encoder/encoder_l2/BiasAdd/ReadVariableOpReadVariableOp9cvae_3_encoder_encoder_l2_biasadd_readvariableop_resource*
_output_shapes
:
*
dtype022
0cvae_3/encoder/encoder_l2/BiasAdd/ReadVariableOp�
!cvae_3/encoder/encoder_l2/BiasAddBiasAdd*cvae_3/encoder/encoder_l2/MatMul:product:08cvae_3/encoder/encoder_l2/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������
2#
!cvae_3/encoder/encoder_l2/BiasAdd�
cvae_3/encoder/encoder_l2/ReluRelu*cvae_3/encoder/encoder_l2/BiasAdd:output:0*
T0*'
_output_shapes
:���������
2 
cvae_3/encoder/encoder_l2/Relu�
/cvae_3/encoder/encoder_l3/MatMul/ReadVariableOpReadVariableOp8cvae_3_encoder_encoder_l3_matmul_readvariableop_resource*
_output_shapes

:
*
dtype021
/cvae_3/encoder/encoder_l3/MatMul/ReadVariableOp�
 cvae_3/encoder/encoder_l3/MatMulMatMul,cvae_3/encoder/encoder_l2/Relu:activations:07cvae_3/encoder/encoder_l3/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2"
 cvae_3/encoder/encoder_l3/MatMul�
0cvae_3/encoder/encoder_l3/BiasAdd/ReadVariableOpReadVariableOp9cvae_3_encoder_encoder_l3_biasadd_readvariableop_resource*
_output_shapes
:*
dtype022
0cvae_3/encoder/encoder_l3/BiasAdd/ReadVariableOp�
!cvae_3/encoder/encoder_l3/BiasAddBiasAdd*cvae_3/encoder/encoder_l3/MatMul:product:08cvae_3/encoder/encoder_l3/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2#
!cvae_3/encoder/encoder_l3/BiasAdd�
cvae_3/encoder/encoder_l3/ReluRelu*cvae_3/encoder/encoder_l3/BiasAdd:output:0*
T0*'
_output_shapes
:���������2 
cvae_3/encoder/encoder_l3/Relu�
/cvae_3/encoder/encoder_l4/MatMul/ReadVariableOpReadVariableOp8cvae_3_encoder_encoder_l4_matmul_readvariableop_resource*
_output_shapes

:*
dtype021
/cvae_3/encoder/encoder_l4/MatMul/ReadVariableOp�
 cvae_3/encoder/encoder_l4/MatMulMatMul,cvae_3/encoder/encoder_l3/Relu:activations:07cvae_3/encoder/encoder_l4/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2"
 cvae_3/encoder/encoder_l4/MatMul�
0cvae_3/encoder/encoder_l4/BiasAdd/ReadVariableOpReadVariableOp9cvae_3_encoder_encoder_l4_biasadd_readvariableop_resource*
_output_shapes
:*
dtype022
0cvae_3/encoder/encoder_l4/BiasAdd/ReadVariableOp�
!cvae_3/encoder/encoder_l4/BiasAddBiasAdd*cvae_3/encoder/encoder_l4/MatMul:product:08cvae_3/encoder/encoder_l4/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2#
!cvae_3/encoder/encoder_l4/BiasAdd�
cvae_3/encoder/encoder_l4/ReluRelu*cvae_3/encoder/encoder_l4/BiasAdd:output:0*
T0*'
_output_shapes
:���������2 
cvae_3/encoder/encoder_l4/Relu�
+cvae_3/encoder/z_mean/MatMul/ReadVariableOpReadVariableOp4cvae_3_encoder_z_mean_matmul_readvariableop_resource*
_output_shapes

:*
dtype02-
+cvae_3/encoder/z_mean/MatMul/ReadVariableOp�
cvae_3/encoder/z_mean/MatMulMatMul,cvae_3/encoder/encoder_l4/Relu:activations:03cvae_3/encoder/z_mean/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2
cvae_3/encoder/z_mean/MatMul�
,cvae_3/encoder/z_mean/BiasAdd/ReadVariableOpReadVariableOp5cvae_3_encoder_z_mean_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02.
,cvae_3/encoder/z_mean/BiasAdd/ReadVariableOp�
cvae_3/encoder/z_mean/BiasAddBiasAdd&cvae_3/encoder/z_mean/MatMul:product:04cvae_3/encoder/z_mean/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2
cvae_3/encoder/z_mean/BiasAdd�
.cvae_3/encoder/z_log_var/MatMul/ReadVariableOpReadVariableOp7cvae_3_encoder_z_log_var_matmul_readvariableop_resource*
_output_shapes

:*
dtype020
.cvae_3/encoder/z_log_var/MatMul/ReadVariableOp�
cvae_3/encoder/z_log_var/MatMulMatMul,cvae_3/encoder/encoder_l4/Relu:activations:06cvae_3/encoder/z_log_var/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2!
cvae_3/encoder/z_log_var/MatMul�
/cvae_3/encoder/z_log_var/BiasAdd/ReadVariableOpReadVariableOp8cvae_3_encoder_z_log_var_biasadd_readvariableop_resource*
_output_shapes
:*
dtype021
/cvae_3/encoder/z_log_var/BiasAdd/ReadVariableOp�
 cvae_3/encoder/z_log_var/BiasAddBiasAdd)cvae_3/encoder/z_log_var/MatMul:product:07cvae_3/encoder/z_log_var/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2"
 cvae_3/encoder/z_log_var/BiasAdd�
cvae_3/encoder/sampling_3/ShapeShape&cvae_3/encoder/z_mean/BiasAdd:output:0*
T0*
_output_shapes
:2!
cvae_3/encoder/sampling_3/Shape�
-cvae_3/encoder/sampling_3/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB: 2/
-cvae_3/encoder/sampling_3/strided_slice/stack�
/cvae_3/encoder/sampling_3/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:21
/cvae_3/encoder/sampling_3/strided_slice/stack_1�
/cvae_3/encoder/sampling_3/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:21
/cvae_3/encoder/sampling_3/strided_slice/stack_2�
'cvae_3/encoder/sampling_3/strided_sliceStridedSlice(cvae_3/encoder/sampling_3/Shape:output:06cvae_3/encoder/sampling_3/strided_slice/stack:output:08cvae_3/encoder/sampling_3/strided_slice/stack_1:output:08cvae_3/encoder/sampling_3/strided_slice/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2)
'cvae_3/encoder/sampling_3/strided_slice�
!cvae_3/encoder/sampling_3/Shape_1Shape&cvae_3/encoder/z_mean/BiasAdd:output:0*
T0*
_output_shapes
:2#
!cvae_3/encoder/sampling_3/Shape_1�
/cvae_3/encoder/sampling_3/strided_slice_1/stackConst*
_output_shapes
:*
dtype0*
valueB:21
/cvae_3/encoder/sampling_3/strided_slice_1/stack�
1cvae_3/encoder/sampling_3/strided_slice_1/stack_1Const*
_output_shapes
:*
dtype0*
valueB:23
1cvae_3/encoder/sampling_3/strided_slice_1/stack_1�
1cvae_3/encoder/sampling_3/strided_slice_1/stack_2Const*
_output_shapes
:*
dtype0*
valueB:23
1cvae_3/encoder/sampling_3/strided_slice_1/stack_2�
)cvae_3/encoder/sampling_3/strided_slice_1StridedSlice*cvae_3/encoder/sampling_3/Shape_1:output:08cvae_3/encoder/sampling_3/strided_slice_1/stack:output:0:cvae_3/encoder/sampling_3/strided_slice_1/stack_1:output:0:cvae_3/encoder/sampling_3/strided_slice_1/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2+
)cvae_3/encoder/sampling_3/strided_slice_1�
-cvae_3/encoder/sampling_3/random_normal/shapePack0cvae_3/encoder/sampling_3/strided_slice:output:02cvae_3/encoder/sampling_3/strided_slice_1:output:0*
N*
T0*
_output_shapes
:2/
-cvae_3/encoder/sampling_3/random_normal/shape�
,cvae_3/encoder/sampling_3/random_normal/meanConst*
_output_shapes
: *
dtype0*
valueB
 *    2.
,cvae_3/encoder/sampling_3/random_normal/mean�
.cvae_3/encoder/sampling_3/random_normal/stddevConst*
_output_shapes
: *
dtype0*
valueB
 *  �?20
.cvae_3/encoder/sampling_3/random_normal/stddev�
<cvae_3/encoder/sampling_3/random_normal/RandomStandardNormalRandomStandardNormal6cvae_3/encoder/sampling_3/random_normal/shape:output:0*
T0*0
_output_shapes
:������������������*
dtype0*
seed���)*
seed2���2>
<cvae_3/encoder/sampling_3/random_normal/RandomStandardNormal�
+cvae_3/encoder/sampling_3/random_normal/mulMulEcvae_3/encoder/sampling_3/random_normal/RandomStandardNormal:output:07cvae_3/encoder/sampling_3/random_normal/stddev:output:0*
T0*0
_output_shapes
:������������������2-
+cvae_3/encoder/sampling_3/random_normal/mul�
'cvae_3/encoder/sampling_3/random_normalAdd/cvae_3/encoder/sampling_3/random_normal/mul:z:05cvae_3/encoder/sampling_3/random_normal/mean:output:0*
T0*0
_output_shapes
:������������������2)
'cvae_3/encoder/sampling_3/random_normal�
cvae_3/encoder/sampling_3/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *   ?2!
cvae_3/encoder/sampling_3/mul/x�
cvae_3/encoder/sampling_3/mulMul(cvae_3/encoder/sampling_3/mul/x:output:0)cvae_3/encoder/z_log_var/BiasAdd:output:0*
T0*'
_output_shapes
:���������2
cvae_3/encoder/sampling_3/mul�
cvae_3/encoder/sampling_3/ExpExp!cvae_3/encoder/sampling_3/mul:z:0*
T0*'
_output_shapes
:���������2
cvae_3/encoder/sampling_3/Exp�
cvae_3/encoder/sampling_3/mul_1Mul!cvae_3/encoder/sampling_3/Exp:y:0+cvae_3/encoder/sampling_3/random_normal:z:0*
T0*'
_output_shapes
:���������2!
cvae_3/encoder/sampling_3/mul_1�
cvae_3/encoder/sampling_3/addAddV2&cvae_3/encoder/z_mean/BiasAdd:output:0#cvae_3/encoder/sampling_3/mul_1:z:0*
T0*'
_output_shapes
:���������2
cvae_3/encoder/sampling_3/add�
cvae_3/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB"    ����2
cvae_3/strided_slice/stack�
cvae_3/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB"        2
cvae_3/strided_slice/stack_1�
cvae_3/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      2
cvae_3/strided_slice/stack_2�
cvae_3/strided_sliceStridedSliceinput_1#cvae_3/strided_slice/stack:output:0%cvae_3/strided_slice/stack_1:output:0%cvae_3/strided_slice/stack_2:output:0*
Index0*
T0*'
_output_shapes
:���������*

begin_mask*
end_mask2
cvae_3/strided_slice�
cvae_3/concatenate/concat/axisConst*
_output_shapes
: *
dtype0*
value	B :2 
cvae_3/concatenate/concat/axis�
cvae_3/concatenate/concatConcatV2!cvae_3/encoder/sampling_3/add:z:0cvae_3/strided_slice:output:0'cvae_3/concatenate/concat/axis:output:0*
N*
T0*'
_output_shapes
:���������2
cvae_3/concatenate/concat�
/cvae_3/decoder/decoder_l2/MatMul/ReadVariableOpReadVariableOp8cvae_3_decoder_decoder_l2_matmul_readvariableop_resource*
_output_shapes

:*
dtype021
/cvae_3/decoder/decoder_l2/MatMul/ReadVariableOp�
 cvae_3/decoder/decoder_l2/MatMulMatMul"cvae_3/concatenate/concat:output:07cvae_3/decoder/decoder_l2/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2"
 cvae_3/decoder/decoder_l2/MatMul�
0cvae_3/decoder/decoder_l2/BiasAdd/ReadVariableOpReadVariableOp9cvae_3_decoder_decoder_l2_biasadd_readvariableop_resource*
_output_shapes
:*
dtype022
0cvae_3/decoder/decoder_l2/BiasAdd/ReadVariableOp�
!cvae_3/decoder/decoder_l2/BiasAddBiasAdd*cvae_3/decoder/decoder_l2/MatMul:product:08cvae_3/decoder/decoder_l2/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2#
!cvae_3/decoder/decoder_l2/BiasAdd�
cvae_3/decoder/decoder_l2/ReluRelu*cvae_3/decoder/decoder_l2/BiasAdd:output:0*
T0*'
_output_shapes
:���������2 
cvae_3/decoder/decoder_l2/Relu�
/cvae_3/decoder/decoder_l3/MatMul/ReadVariableOpReadVariableOp8cvae_3_decoder_decoder_l3_matmul_readvariableop_resource*
_output_shapes

:*
dtype021
/cvae_3/decoder/decoder_l3/MatMul/ReadVariableOp�
 cvae_3/decoder/decoder_l3/MatMulMatMul,cvae_3/decoder/decoder_l2/Relu:activations:07cvae_3/decoder/decoder_l3/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2"
 cvae_3/decoder/decoder_l3/MatMul�
0cvae_3/decoder/decoder_l3/BiasAdd/ReadVariableOpReadVariableOp9cvae_3_decoder_decoder_l3_biasadd_readvariableop_resource*
_output_shapes
:*
dtype022
0cvae_3/decoder/decoder_l3/BiasAdd/ReadVariableOp�
!cvae_3/decoder/decoder_l3/BiasAddBiasAdd*cvae_3/decoder/decoder_l3/MatMul:product:08cvae_3/decoder/decoder_l3/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2#
!cvae_3/decoder/decoder_l3/BiasAdd�
cvae_3/decoder/decoder_l3/ReluRelu*cvae_3/decoder/decoder_l3/BiasAdd:output:0*
T0*'
_output_shapes
:���������2 
cvae_3/decoder/decoder_l3/Relu�
/cvae_3/decoder/decoder_l4/MatMul/ReadVariableOpReadVariableOp8cvae_3_decoder_decoder_l4_matmul_readvariableop_resource*
_output_shapes

:
*
dtype021
/cvae_3/decoder/decoder_l4/MatMul/ReadVariableOp�
 cvae_3/decoder/decoder_l4/MatMulMatMul,cvae_3/decoder/decoder_l3/Relu:activations:07cvae_3/decoder/decoder_l4/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������
2"
 cvae_3/decoder/decoder_l4/MatMul�
0cvae_3/decoder/decoder_l4/BiasAdd/ReadVariableOpReadVariableOp9cvae_3_decoder_decoder_l4_biasadd_readvariableop_resource*
_output_shapes
:
*
dtype022
0cvae_3/decoder/decoder_l4/BiasAdd/ReadVariableOp�
!cvae_3/decoder/decoder_l4/BiasAddBiasAdd*cvae_3/decoder/decoder_l4/MatMul:product:08cvae_3/decoder/decoder_l4/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������
2#
!cvae_3/decoder/decoder_l4/BiasAdd�
cvae_3/decoder/decoder_l4/ReluRelu*cvae_3/decoder/decoder_l4/BiasAdd:output:0*
T0*'
_output_shapes
:���������
2 
cvae_3/decoder/decoder_l4/Relu�
1cvae_3/decoder/output_layer/MatMul/ReadVariableOpReadVariableOp:cvae_3_decoder_output_layer_matmul_readvariableop_resource*
_output_shapes

:
*
dtype023
1cvae_3/decoder/output_layer/MatMul/ReadVariableOp�
"cvae_3/decoder/output_layer/MatMulMatMul,cvae_3/decoder/decoder_l4/Relu:activations:09cvae_3/decoder/output_layer/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2$
"cvae_3/decoder/output_layer/MatMul�
2cvae_3/decoder/output_layer/BiasAdd/ReadVariableOpReadVariableOp;cvae_3_decoder_output_layer_biasadd_readvariableop_resource*
_output_shapes
:*
dtype024
2cvae_3/decoder/output_layer/BiasAdd/ReadVariableOp�
#cvae_3/decoder/output_layer/BiasAddBiasAdd,cvae_3/decoder/output_layer/MatMul:product:0:cvae_3/decoder/output_layer/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2%
#cvae_3/decoder/output_layer/BiasAdd�
 cvae_3/decoder/output_layer/ReluRelu,cvae_3/decoder/output_layer/BiasAdd:output:0*
T0*'
_output_shapes
:���������2"
 cvae_3/decoder/output_layer/Relu�
IdentityIdentity.cvae_3/decoder/output_layer/Relu:activations:01^cvae_3/decoder/decoder_l2/BiasAdd/ReadVariableOp0^cvae_3/decoder/decoder_l2/MatMul/ReadVariableOp1^cvae_3/decoder/decoder_l3/BiasAdd/ReadVariableOp0^cvae_3/decoder/decoder_l3/MatMul/ReadVariableOp1^cvae_3/decoder/decoder_l4/BiasAdd/ReadVariableOp0^cvae_3/decoder/decoder_l4/MatMul/ReadVariableOp3^cvae_3/decoder/output_layer/BiasAdd/ReadVariableOp2^cvae_3/decoder/output_layer/MatMul/ReadVariableOp1^cvae_3/encoder/encoder_l2/BiasAdd/ReadVariableOp0^cvae_3/encoder/encoder_l2/MatMul/ReadVariableOp1^cvae_3/encoder/encoder_l3/BiasAdd/ReadVariableOp0^cvae_3/encoder/encoder_l3/MatMul/ReadVariableOp1^cvae_3/encoder/encoder_l4/BiasAdd/ReadVariableOp0^cvae_3/encoder/encoder_l4/MatMul/ReadVariableOp0^cvae_3/encoder/z_log_var/BiasAdd/ReadVariableOp/^cvae_3/encoder/z_log_var/MatMul/ReadVariableOp-^cvae_3/encoder/z_mean/BiasAdd/ReadVariableOp,^cvae_3/encoder/z_mean/MatMul/ReadVariableOp*
T0*'
_output_shapes
:���������2

Identity"
identityIdentity:output:0*n
_input_shapes]
[:���������::::::::::::::::::2d
0cvae_3/decoder/decoder_l2/BiasAdd/ReadVariableOp0cvae_3/decoder/decoder_l2/BiasAdd/ReadVariableOp2b
/cvae_3/decoder/decoder_l2/MatMul/ReadVariableOp/cvae_3/decoder/decoder_l2/MatMul/ReadVariableOp2d
0cvae_3/decoder/decoder_l3/BiasAdd/ReadVariableOp0cvae_3/decoder/decoder_l3/BiasAdd/ReadVariableOp2b
/cvae_3/decoder/decoder_l3/MatMul/ReadVariableOp/cvae_3/decoder/decoder_l3/MatMul/ReadVariableOp2d
0cvae_3/decoder/decoder_l4/BiasAdd/ReadVariableOp0cvae_3/decoder/decoder_l4/BiasAdd/ReadVariableOp2b
/cvae_3/decoder/decoder_l4/MatMul/ReadVariableOp/cvae_3/decoder/decoder_l4/MatMul/ReadVariableOp2h
2cvae_3/decoder/output_layer/BiasAdd/ReadVariableOp2cvae_3/decoder/output_layer/BiasAdd/ReadVariableOp2f
1cvae_3/decoder/output_layer/MatMul/ReadVariableOp1cvae_3/decoder/output_layer/MatMul/ReadVariableOp2d
0cvae_3/encoder/encoder_l2/BiasAdd/ReadVariableOp0cvae_3/encoder/encoder_l2/BiasAdd/ReadVariableOp2b
/cvae_3/encoder/encoder_l2/MatMul/ReadVariableOp/cvae_3/encoder/encoder_l2/MatMul/ReadVariableOp2d
0cvae_3/encoder/encoder_l3/BiasAdd/ReadVariableOp0cvae_3/encoder/encoder_l3/BiasAdd/ReadVariableOp2b
/cvae_3/encoder/encoder_l3/MatMul/ReadVariableOp/cvae_3/encoder/encoder_l3/MatMul/ReadVariableOp2d
0cvae_3/encoder/encoder_l4/BiasAdd/ReadVariableOp0cvae_3/encoder/encoder_l4/BiasAdd/ReadVariableOp2b
/cvae_3/encoder/encoder_l4/MatMul/ReadVariableOp/cvae_3/encoder/encoder_l4/MatMul/ReadVariableOp2b
/cvae_3/encoder/z_log_var/BiasAdd/ReadVariableOp/cvae_3/encoder/z_log_var/BiasAdd/ReadVariableOp2`
.cvae_3/encoder/z_log_var/MatMul/ReadVariableOp.cvae_3/encoder/z_log_var/MatMul/ReadVariableOp2\
,cvae_3/encoder/z_mean/BiasAdd/ReadVariableOp,cvae_3/encoder/z_mean/BiasAdd/ReadVariableOp2Z
+cvae_3/encoder/z_mean/MatMul/ReadVariableOp+cvae_3/encoder/z_mean/MatMul/ReadVariableOp:P L
'
_output_shapes
:���������
!
_user_specified_name	input_1
�
�
+__inference_encoder_l2_layer_call_fn_142931

inputs
unknown
	unknown_0
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������
*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *O
fJRH
F__inference_encoder_l2_layer_call_and_return_conditional_losses_1411672
StatefulPartitionedCall�
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*'
_output_shapes
:���������
2

Identity"
identityIdentity:output:0*.
_input_shapes
:���������::22
StatefulPartitionedCallStatefulPartitionedCall:O K
'
_output_shapes
:���������
 
_user_specified_nameinputs
�
�
(__inference_encoder_layer_call_fn_142805

inputs
unknown
	unknown_0
	unknown_1
	unknown_2
	unknown_3
	unknown_4
	unknown_5
	unknown_6
	unknown_7
	unknown_8
identity

identity_1

identity_2��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6	unknown_7	unknown_8*
Tin
2*
Tout
2*
_collective_manager_ids
 *M
_output_shapes;
9:���������:���������:���������*,
_read_only_resource_inputs

	
*-
config_proto

CPU

GPU 2J 8� *L
fGRE
C__inference_encoder_layer_call_and_return_conditional_losses_1414552
StatefulPartitionedCall�
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*'
_output_shapes
:���������2

Identity�

Identity_1Identity StatefulPartitionedCall:output:1^StatefulPartitionedCall*
T0*'
_output_shapes
:���������2

Identity_1�

Identity_2Identity StatefulPartitionedCall:output:2^StatefulPartitionedCall*
T0*'
_output_shapes
:���������2

Identity_2"
identityIdentity:output:0"!

identity_1Identity_1:output:0"!

identity_2Identity_2:output:0*N
_input_shapes=
;:���������::::::::::22
StatefulPartitionedCallStatefulPartitionedCall:O K
'
_output_shapes
:���������
 
_user_specified_nameinputs
��
� 
"__inference__traced_restore_143520
file_prefix
assignvariableop_adam_iter"
assignvariableop_1_adam_beta_1"
assignvariableop_2_adam_beta_2!
assignvariableop_3_adam_decay)
%assignvariableop_4_adam_learning_rate(
$assignvariableop_5_encoder_l2_kernel&
"assignvariableop_6_encoder_l2_bias(
$assignvariableop_7_encoder_l3_kernel&
"assignvariableop_8_encoder_l3_bias(
$assignvariableop_9_encoder_l4_kernel'
#assignvariableop_10_encoder_l4_bias%
!assignvariableop_11_z_mean_kernel#
assignvariableop_12_z_mean_bias(
$assignvariableop_13_z_log_var_kernel&
"assignvariableop_14_z_log_var_bias)
%assignvariableop_15_decoder_l2_kernel'
#assignvariableop_16_decoder_l2_bias)
%assignvariableop_17_decoder_l3_kernel'
#assignvariableop_18_decoder_l3_bias)
%assignvariableop_19_decoder_l4_kernel'
#assignvariableop_20_decoder_l4_bias+
'assignvariableop_21_output_layer_kernel)
%assignvariableop_22_output_layer_bias
assignvariableop_23_total
assignvariableop_24_count0
,assignvariableop_25_adam_encoder_l2_kernel_m.
*assignvariableop_26_adam_encoder_l2_bias_m0
,assignvariableop_27_adam_encoder_l3_kernel_m.
*assignvariableop_28_adam_encoder_l3_bias_m0
,assignvariableop_29_adam_encoder_l4_kernel_m.
*assignvariableop_30_adam_encoder_l4_bias_m,
(assignvariableop_31_adam_z_mean_kernel_m*
&assignvariableop_32_adam_z_mean_bias_m/
+assignvariableop_33_adam_z_log_var_kernel_m-
)assignvariableop_34_adam_z_log_var_bias_m0
,assignvariableop_35_adam_decoder_l2_kernel_m.
*assignvariableop_36_adam_decoder_l2_bias_m0
,assignvariableop_37_adam_decoder_l3_kernel_m.
*assignvariableop_38_adam_decoder_l3_bias_m0
,assignvariableop_39_adam_decoder_l4_kernel_m.
*assignvariableop_40_adam_decoder_l4_bias_m2
.assignvariableop_41_adam_output_layer_kernel_m0
,assignvariableop_42_adam_output_layer_bias_m0
,assignvariableop_43_adam_encoder_l2_kernel_v.
*assignvariableop_44_adam_encoder_l2_bias_v0
,assignvariableop_45_adam_encoder_l3_kernel_v.
*assignvariableop_46_adam_encoder_l3_bias_v0
,assignvariableop_47_adam_encoder_l4_kernel_v.
*assignvariableop_48_adam_encoder_l4_bias_v,
(assignvariableop_49_adam_z_mean_kernel_v*
&assignvariableop_50_adam_z_mean_bias_v/
+assignvariableop_51_adam_z_log_var_kernel_v-
)assignvariableop_52_adam_z_log_var_bias_v0
,assignvariableop_53_adam_decoder_l2_kernel_v.
*assignvariableop_54_adam_decoder_l2_bias_v0
,assignvariableop_55_adam_decoder_l3_kernel_v.
*assignvariableop_56_adam_decoder_l3_bias_v0
,assignvariableop_57_adam_decoder_l4_kernel_v.
*assignvariableop_58_adam_decoder_l4_bias_v2
.assignvariableop_59_adam_output_layer_kernel_v0
,assignvariableop_60_adam_output_layer_bias_v
identity_62��AssignVariableOp�AssignVariableOp_1�AssignVariableOp_10�AssignVariableOp_11�AssignVariableOp_12�AssignVariableOp_13�AssignVariableOp_14�AssignVariableOp_15�AssignVariableOp_16�AssignVariableOp_17�AssignVariableOp_18�AssignVariableOp_19�AssignVariableOp_2�AssignVariableOp_20�AssignVariableOp_21�AssignVariableOp_22�AssignVariableOp_23�AssignVariableOp_24�AssignVariableOp_25�AssignVariableOp_26�AssignVariableOp_27�AssignVariableOp_28�AssignVariableOp_29�AssignVariableOp_3�AssignVariableOp_30�AssignVariableOp_31�AssignVariableOp_32�AssignVariableOp_33�AssignVariableOp_34�AssignVariableOp_35�AssignVariableOp_36�AssignVariableOp_37�AssignVariableOp_38�AssignVariableOp_39�AssignVariableOp_4�AssignVariableOp_40�AssignVariableOp_41�AssignVariableOp_42�AssignVariableOp_43�AssignVariableOp_44�AssignVariableOp_45�AssignVariableOp_46�AssignVariableOp_47�AssignVariableOp_48�AssignVariableOp_49�AssignVariableOp_5�AssignVariableOp_50�AssignVariableOp_51�AssignVariableOp_52�AssignVariableOp_53�AssignVariableOp_54�AssignVariableOp_55�AssignVariableOp_56�AssignVariableOp_57�AssignVariableOp_58�AssignVariableOp_59�AssignVariableOp_6�AssignVariableOp_60�AssignVariableOp_7�AssignVariableOp_8�AssignVariableOp_9�
RestoreV2/tensor_namesConst"/device:CPU:0*
_output_shapes
:>*
dtype0*�
value�B�>B)optimizer/iter/.ATTRIBUTES/VARIABLE_VALUEB+optimizer/beta_1/.ATTRIBUTES/VARIABLE_VALUEB+optimizer/beta_2/.ATTRIBUTES/VARIABLE_VALUEB*optimizer/decay/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/learning_rate/.ATTRIBUTES/VARIABLE_VALUEB&variables/0/.ATTRIBUTES/VARIABLE_VALUEB&variables/1/.ATTRIBUTES/VARIABLE_VALUEB&variables/2/.ATTRIBUTES/VARIABLE_VALUEB&variables/3/.ATTRIBUTES/VARIABLE_VALUEB&variables/4/.ATTRIBUTES/VARIABLE_VALUEB&variables/5/.ATTRIBUTES/VARIABLE_VALUEB&variables/6/.ATTRIBUTES/VARIABLE_VALUEB&variables/7/.ATTRIBUTES/VARIABLE_VALUEB&variables/8/.ATTRIBUTES/VARIABLE_VALUEB&variables/9/.ATTRIBUTES/VARIABLE_VALUEB'variables/10/.ATTRIBUTES/VARIABLE_VALUEB'variables/11/.ATTRIBUTES/VARIABLE_VALUEB'variables/12/.ATTRIBUTES/VARIABLE_VALUEB'variables/13/.ATTRIBUTES/VARIABLE_VALUEB'variables/14/.ATTRIBUTES/VARIABLE_VALUEB'variables/15/.ATTRIBUTES/VARIABLE_VALUEB'variables/16/.ATTRIBUTES/VARIABLE_VALUEB'variables/17/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/0/total/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/0/count/.ATTRIBUTES/VARIABLE_VALUEBBvariables/0/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBBvariables/1/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBBvariables/2/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBBvariables/3/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBBvariables/4/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBBvariables/5/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBBvariables/6/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBBvariables/7/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBBvariables/8/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBBvariables/9/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBCvariables/10/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBCvariables/11/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBCvariables/12/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBCvariables/13/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBCvariables/14/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBCvariables/15/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBCvariables/16/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBCvariables/17/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBBvariables/0/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBBvariables/1/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBBvariables/2/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBBvariables/3/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBBvariables/4/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBBvariables/5/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBBvariables/6/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBBvariables/7/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBBvariables/8/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBBvariables/9/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBCvariables/10/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBCvariables/11/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBCvariables/12/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBCvariables/13/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBCvariables/14/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBCvariables/15/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBCvariables/16/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBCvariables/17/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEB_CHECKPOINTABLE_OBJECT_GRAPH2
RestoreV2/tensor_names�
RestoreV2/shape_and_slicesConst"/device:CPU:0*
_output_shapes
:>*
dtype0*�
value�B�>B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B 2
RestoreV2/shape_and_slices�
	RestoreV2	RestoreV2file_prefixRestoreV2/tensor_names:output:0#RestoreV2/shape_and_slices:output:0"/device:CPU:0*�
_output_shapes�
�::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*L
dtypesB
@2>	2
	RestoreV2g
IdentityIdentityRestoreV2:tensors:0"/device:CPU:0*
T0	*
_output_shapes
:2

Identity�
AssignVariableOpAssignVariableOpassignvariableop_adam_iterIdentity:output:0"/device:CPU:0*
_output_shapes
 *
dtype0	2
AssignVariableOpk

Identity_1IdentityRestoreV2:tensors:1"/device:CPU:0*
T0*
_output_shapes
:2

Identity_1�
AssignVariableOp_1AssignVariableOpassignvariableop_1_adam_beta_1Identity_1:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_1k

Identity_2IdentityRestoreV2:tensors:2"/device:CPU:0*
T0*
_output_shapes
:2

Identity_2�
AssignVariableOp_2AssignVariableOpassignvariableop_2_adam_beta_2Identity_2:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_2k

Identity_3IdentityRestoreV2:tensors:3"/device:CPU:0*
T0*
_output_shapes
:2

Identity_3�
AssignVariableOp_3AssignVariableOpassignvariableop_3_adam_decayIdentity_3:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_3k

Identity_4IdentityRestoreV2:tensors:4"/device:CPU:0*
T0*
_output_shapes
:2

Identity_4�
AssignVariableOp_4AssignVariableOp%assignvariableop_4_adam_learning_rateIdentity_4:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_4k

Identity_5IdentityRestoreV2:tensors:5"/device:CPU:0*
T0*
_output_shapes
:2

Identity_5�
AssignVariableOp_5AssignVariableOp$assignvariableop_5_encoder_l2_kernelIdentity_5:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_5k

Identity_6IdentityRestoreV2:tensors:6"/device:CPU:0*
T0*
_output_shapes
:2

Identity_6�
AssignVariableOp_6AssignVariableOp"assignvariableop_6_encoder_l2_biasIdentity_6:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_6k

Identity_7IdentityRestoreV2:tensors:7"/device:CPU:0*
T0*
_output_shapes
:2

Identity_7�
AssignVariableOp_7AssignVariableOp$assignvariableop_7_encoder_l3_kernelIdentity_7:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_7k

Identity_8IdentityRestoreV2:tensors:8"/device:CPU:0*
T0*
_output_shapes
:2

Identity_8�
AssignVariableOp_8AssignVariableOp"assignvariableop_8_encoder_l3_biasIdentity_8:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_8k

Identity_9IdentityRestoreV2:tensors:9"/device:CPU:0*
T0*
_output_shapes
:2

Identity_9�
AssignVariableOp_9AssignVariableOp$assignvariableop_9_encoder_l4_kernelIdentity_9:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_9n
Identity_10IdentityRestoreV2:tensors:10"/device:CPU:0*
T0*
_output_shapes
:2
Identity_10�
AssignVariableOp_10AssignVariableOp#assignvariableop_10_encoder_l4_biasIdentity_10:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_10n
Identity_11IdentityRestoreV2:tensors:11"/device:CPU:0*
T0*
_output_shapes
:2
Identity_11�
AssignVariableOp_11AssignVariableOp!assignvariableop_11_z_mean_kernelIdentity_11:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_11n
Identity_12IdentityRestoreV2:tensors:12"/device:CPU:0*
T0*
_output_shapes
:2
Identity_12�
AssignVariableOp_12AssignVariableOpassignvariableop_12_z_mean_biasIdentity_12:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_12n
Identity_13IdentityRestoreV2:tensors:13"/device:CPU:0*
T0*
_output_shapes
:2
Identity_13�
AssignVariableOp_13AssignVariableOp$assignvariableop_13_z_log_var_kernelIdentity_13:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_13n
Identity_14IdentityRestoreV2:tensors:14"/device:CPU:0*
T0*
_output_shapes
:2
Identity_14�
AssignVariableOp_14AssignVariableOp"assignvariableop_14_z_log_var_biasIdentity_14:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_14n
Identity_15IdentityRestoreV2:tensors:15"/device:CPU:0*
T0*
_output_shapes
:2
Identity_15�
AssignVariableOp_15AssignVariableOp%assignvariableop_15_decoder_l2_kernelIdentity_15:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_15n
Identity_16IdentityRestoreV2:tensors:16"/device:CPU:0*
T0*
_output_shapes
:2
Identity_16�
AssignVariableOp_16AssignVariableOp#assignvariableop_16_decoder_l2_biasIdentity_16:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_16n
Identity_17IdentityRestoreV2:tensors:17"/device:CPU:0*
T0*
_output_shapes
:2
Identity_17�
AssignVariableOp_17AssignVariableOp%assignvariableop_17_decoder_l3_kernelIdentity_17:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_17n
Identity_18IdentityRestoreV2:tensors:18"/device:CPU:0*
T0*
_output_shapes
:2
Identity_18�
AssignVariableOp_18AssignVariableOp#assignvariableop_18_decoder_l3_biasIdentity_18:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_18n
Identity_19IdentityRestoreV2:tensors:19"/device:CPU:0*
T0*
_output_shapes
:2
Identity_19�
AssignVariableOp_19AssignVariableOp%assignvariableop_19_decoder_l4_kernelIdentity_19:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_19n
Identity_20IdentityRestoreV2:tensors:20"/device:CPU:0*
T0*
_output_shapes
:2
Identity_20�
AssignVariableOp_20AssignVariableOp#assignvariableop_20_decoder_l4_biasIdentity_20:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_20n
Identity_21IdentityRestoreV2:tensors:21"/device:CPU:0*
T0*
_output_shapes
:2
Identity_21�
AssignVariableOp_21AssignVariableOp'assignvariableop_21_output_layer_kernelIdentity_21:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_21n
Identity_22IdentityRestoreV2:tensors:22"/device:CPU:0*
T0*
_output_shapes
:2
Identity_22�
AssignVariableOp_22AssignVariableOp%assignvariableop_22_output_layer_biasIdentity_22:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_22n
Identity_23IdentityRestoreV2:tensors:23"/device:CPU:0*
T0*
_output_shapes
:2
Identity_23�
AssignVariableOp_23AssignVariableOpassignvariableop_23_totalIdentity_23:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_23n
Identity_24IdentityRestoreV2:tensors:24"/device:CPU:0*
T0*
_output_shapes
:2
Identity_24�
AssignVariableOp_24AssignVariableOpassignvariableop_24_countIdentity_24:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_24n
Identity_25IdentityRestoreV2:tensors:25"/device:CPU:0*
T0*
_output_shapes
:2
Identity_25�
AssignVariableOp_25AssignVariableOp,assignvariableop_25_adam_encoder_l2_kernel_mIdentity_25:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_25n
Identity_26IdentityRestoreV2:tensors:26"/device:CPU:0*
T0*
_output_shapes
:2
Identity_26�
AssignVariableOp_26AssignVariableOp*assignvariableop_26_adam_encoder_l2_bias_mIdentity_26:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_26n
Identity_27IdentityRestoreV2:tensors:27"/device:CPU:0*
T0*
_output_shapes
:2
Identity_27�
AssignVariableOp_27AssignVariableOp,assignvariableop_27_adam_encoder_l3_kernel_mIdentity_27:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_27n
Identity_28IdentityRestoreV2:tensors:28"/device:CPU:0*
T0*
_output_shapes
:2
Identity_28�
AssignVariableOp_28AssignVariableOp*assignvariableop_28_adam_encoder_l3_bias_mIdentity_28:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_28n
Identity_29IdentityRestoreV2:tensors:29"/device:CPU:0*
T0*
_output_shapes
:2
Identity_29�
AssignVariableOp_29AssignVariableOp,assignvariableop_29_adam_encoder_l4_kernel_mIdentity_29:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_29n
Identity_30IdentityRestoreV2:tensors:30"/device:CPU:0*
T0*
_output_shapes
:2
Identity_30�
AssignVariableOp_30AssignVariableOp*assignvariableop_30_adam_encoder_l4_bias_mIdentity_30:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_30n
Identity_31IdentityRestoreV2:tensors:31"/device:CPU:0*
T0*
_output_shapes
:2
Identity_31�
AssignVariableOp_31AssignVariableOp(assignvariableop_31_adam_z_mean_kernel_mIdentity_31:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_31n
Identity_32IdentityRestoreV2:tensors:32"/device:CPU:0*
T0*
_output_shapes
:2
Identity_32�
AssignVariableOp_32AssignVariableOp&assignvariableop_32_adam_z_mean_bias_mIdentity_32:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_32n
Identity_33IdentityRestoreV2:tensors:33"/device:CPU:0*
T0*
_output_shapes
:2
Identity_33�
AssignVariableOp_33AssignVariableOp+assignvariableop_33_adam_z_log_var_kernel_mIdentity_33:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_33n
Identity_34IdentityRestoreV2:tensors:34"/device:CPU:0*
T0*
_output_shapes
:2
Identity_34�
AssignVariableOp_34AssignVariableOp)assignvariableop_34_adam_z_log_var_bias_mIdentity_34:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_34n
Identity_35IdentityRestoreV2:tensors:35"/device:CPU:0*
T0*
_output_shapes
:2
Identity_35�
AssignVariableOp_35AssignVariableOp,assignvariableop_35_adam_decoder_l2_kernel_mIdentity_35:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_35n
Identity_36IdentityRestoreV2:tensors:36"/device:CPU:0*
T0*
_output_shapes
:2
Identity_36�
AssignVariableOp_36AssignVariableOp*assignvariableop_36_adam_decoder_l2_bias_mIdentity_36:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_36n
Identity_37IdentityRestoreV2:tensors:37"/device:CPU:0*
T0*
_output_shapes
:2
Identity_37�
AssignVariableOp_37AssignVariableOp,assignvariableop_37_adam_decoder_l3_kernel_mIdentity_37:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_37n
Identity_38IdentityRestoreV2:tensors:38"/device:CPU:0*
T0*
_output_shapes
:2
Identity_38�
AssignVariableOp_38AssignVariableOp*assignvariableop_38_adam_decoder_l3_bias_mIdentity_38:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_38n
Identity_39IdentityRestoreV2:tensors:39"/device:CPU:0*
T0*
_output_shapes
:2
Identity_39�
AssignVariableOp_39AssignVariableOp,assignvariableop_39_adam_decoder_l4_kernel_mIdentity_39:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_39n
Identity_40IdentityRestoreV2:tensors:40"/device:CPU:0*
T0*
_output_shapes
:2
Identity_40�
AssignVariableOp_40AssignVariableOp*assignvariableop_40_adam_decoder_l4_bias_mIdentity_40:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_40n
Identity_41IdentityRestoreV2:tensors:41"/device:CPU:0*
T0*
_output_shapes
:2
Identity_41�
AssignVariableOp_41AssignVariableOp.assignvariableop_41_adam_output_layer_kernel_mIdentity_41:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_41n
Identity_42IdentityRestoreV2:tensors:42"/device:CPU:0*
T0*
_output_shapes
:2
Identity_42�
AssignVariableOp_42AssignVariableOp,assignvariableop_42_adam_output_layer_bias_mIdentity_42:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_42n
Identity_43IdentityRestoreV2:tensors:43"/device:CPU:0*
T0*
_output_shapes
:2
Identity_43�
AssignVariableOp_43AssignVariableOp,assignvariableop_43_adam_encoder_l2_kernel_vIdentity_43:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_43n
Identity_44IdentityRestoreV2:tensors:44"/device:CPU:0*
T0*
_output_shapes
:2
Identity_44�
AssignVariableOp_44AssignVariableOp*assignvariableop_44_adam_encoder_l2_bias_vIdentity_44:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_44n
Identity_45IdentityRestoreV2:tensors:45"/device:CPU:0*
T0*
_output_shapes
:2
Identity_45�
AssignVariableOp_45AssignVariableOp,assignvariableop_45_adam_encoder_l3_kernel_vIdentity_45:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_45n
Identity_46IdentityRestoreV2:tensors:46"/device:CPU:0*
T0*
_output_shapes
:2
Identity_46�
AssignVariableOp_46AssignVariableOp*assignvariableop_46_adam_encoder_l3_bias_vIdentity_46:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_46n
Identity_47IdentityRestoreV2:tensors:47"/device:CPU:0*
T0*
_output_shapes
:2
Identity_47�
AssignVariableOp_47AssignVariableOp,assignvariableop_47_adam_encoder_l4_kernel_vIdentity_47:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_47n
Identity_48IdentityRestoreV2:tensors:48"/device:CPU:0*
T0*
_output_shapes
:2
Identity_48�
AssignVariableOp_48AssignVariableOp*assignvariableop_48_adam_encoder_l4_bias_vIdentity_48:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_48n
Identity_49IdentityRestoreV2:tensors:49"/device:CPU:0*
T0*
_output_shapes
:2
Identity_49�
AssignVariableOp_49AssignVariableOp(assignvariableop_49_adam_z_mean_kernel_vIdentity_49:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_49n
Identity_50IdentityRestoreV2:tensors:50"/device:CPU:0*
T0*
_output_shapes
:2
Identity_50�
AssignVariableOp_50AssignVariableOp&assignvariableop_50_adam_z_mean_bias_vIdentity_50:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_50n
Identity_51IdentityRestoreV2:tensors:51"/device:CPU:0*
T0*
_output_shapes
:2
Identity_51�
AssignVariableOp_51AssignVariableOp+assignvariableop_51_adam_z_log_var_kernel_vIdentity_51:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_51n
Identity_52IdentityRestoreV2:tensors:52"/device:CPU:0*
T0*
_output_shapes
:2
Identity_52�
AssignVariableOp_52AssignVariableOp)assignvariableop_52_adam_z_log_var_bias_vIdentity_52:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_52n
Identity_53IdentityRestoreV2:tensors:53"/device:CPU:0*
T0*
_output_shapes
:2
Identity_53�
AssignVariableOp_53AssignVariableOp,assignvariableop_53_adam_decoder_l2_kernel_vIdentity_53:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_53n
Identity_54IdentityRestoreV2:tensors:54"/device:CPU:0*
T0*
_output_shapes
:2
Identity_54�
AssignVariableOp_54AssignVariableOp*assignvariableop_54_adam_decoder_l2_bias_vIdentity_54:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_54n
Identity_55IdentityRestoreV2:tensors:55"/device:CPU:0*
T0*
_output_shapes
:2
Identity_55�
AssignVariableOp_55AssignVariableOp,assignvariableop_55_adam_decoder_l3_kernel_vIdentity_55:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_55n
Identity_56IdentityRestoreV2:tensors:56"/device:CPU:0*
T0*
_output_shapes
:2
Identity_56�
AssignVariableOp_56AssignVariableOp*assignvariableop_56_adam_decoder_l3_bias_vIdentity_56:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_56n
Identity_57IdentityRestoreV2:tensors:57"/device:CPU:0*
T0*
_output_shapes
:2
Identity_57�
AssignVariableOp_57AssignVariableOp,assignvariableop_57_adam_decoder_l4_kernel_vIdentity_57:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_57n
Identity_58IdentityRestoreV2:tensors:58"/device:CPU:0*
T0*
_output_shapes
:2
Identity_58�
AssignVariableOp_58AssignVariableOp*assignvariableop_58_adam_decoder_l4_bias_vIdentity_58:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_58n
Identity_59IdentityRestoreV2:tensors:59"/device:CPU:0*
T0*
_output_shapes
:2
Identity_59�
AssignVariableOp_59AssignVariableOp.assignvariableop_59_adam_output_layer_kernel_vIdentity_59:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_59n
Identity_60IdentityRestoreV2:tensors:60"/device:CPU:0*
T0*
_output_shapes
:2
Identity_60�
AssignVariableOp_60AssignVariableOp,assignvariableop_60_adam_output_layer_bias_vIdentity_60:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_609
NoOpNoOp"/device:CPU:0*
_output_shapes
 2
NoOp�
Identity_61Identityfile_prefix^AssignVariableOp^AssignVariableOp_1^AssignVariableOp_10^AssignVariableOp_11^AssignVariableOp_12^AssignVariableOp_13^AssignVariableOp_14^AssignVariableOp_15^AssignVariableOp_16^AssignVariableOp_17^AssignVariableOp_18^AssignVariableOp_19^AssignVariableOp_2^AssignVariableOp_20^AssignVariableOp_21^AssignVariableOp_22^AssignVariableOp_23^AssignVariableOp_24^AssignVariableOp_25^AssignVariableOp_26^AssignVariableOp_27^AssignVariableOp_28^AssignVariableOp_29^AssignVariableOp_3^AssignVariableOp_30^AssignVariableOp_31^AssignVariableOp_32^AssignVariableOp_33^AssignVariableOp_34^AssignVariableOp_35^AssignVariableOp_36^AssignVariableOp_37^AssignVariableOp_38^AssignVariableOp_39^AssignVariableOp_4^AssignVariableOp_40^AssignVariableOp_41^AssignVariableOp_42^AssignVariableOp_43^AssignVariableOp_44^AssignVariableOp_45^AssignVariableOp_46^AssignVariableOp_47^AssignVariableOp_48^AssignVariableOp_49^AssignVariableOp_5^AssignVariableOp_50^AssignVariableOp_51^AssignVariableOp_52^AssignVariableOp_53^AssignVariableOp_54^AssignVariableOp_55^AssignVariableOp_56^AssignVariableOp_57^AssignVariableOp_58^AssignVariableOp_59^AssignVariableOp_6^AssignVariableOp_60^AssignVariableOp_7^AssignVariableOp_8^AssignVariableOp_9^NoOp"/device:CPU:0*
T0*
_output_shapes
: 2
Identity_61�
Identity_62IdentityIdentity_61:output:0^AssignVariableOp^AssignVariableOp_1^AssignVariableOp_10^AssignVariableOp_11^AssignVariableOp_12^AssignVariableOp_13^AssignVariableOp_14^AssignVariableOp_15^AssignVariableOp_16^AssignVariableOp_17^AssignVariableOp_18^AssignVariableOp_19^AssignVariableOp_2^AssignVariableOp_20^AssignVariableOp_21^AssignVariableOp_22^AssignVariableOp_23^AssignVariableOp_24^AssignVariableOp_25^AssignVariableOp_26^AssignVariableOp_27^AssignVariableOp_28^AssignVariableOp_29^AssignVariableOp_3^AssignVariableOp_30^AssignVariableOp_31^AssignVariableOp_32^AssignVariableOp_33^AssignVariableOp_34^AssignVariableOp_35^AssignVariableOp_36^AssignVariableOp_37^AssignVariableOp_38^AssignVariableOp_39^AssignVariableOp_4^AssignVariableOp_40^AssignVariableOp_41^AssignVariableOp_42^AssignVariableOp_43^AssignVariableOp_44^AssignVariableOp_45^AssignVariableOp_46^AssignVariableOp_47^AssignVariableOp_48^AssignVariableOp_49^AssignVariableOp_5^AssignVariableOp_50^AssignVariableOp_51^AssignVariableOp_52^AssignVariableOp_53^AssignVariableOp_54^AssignVariableOp_55^AssignVariableOp_56^AssignVariableOp_57^AssignVariableOp_58^AssignVariableOp_59^AssignVariableOp_6^AssignVariableOp_60^AssignVariableOp_7^AssignVariableOp_8^AssignVariableOp_9*
T0*
_output_shapes
: 2
Identity_62"#
identity_62Identity_62:output:0*�
_input_shapes�
�: :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::2$
AssignVariableOpAssignVariableOp2(
AssignVariableOp_1AssignVariableOp_12*
AssignVariableOp_10AssignVariableOp_102*
AssignVariableOp_11AssignVariableOp_112*
AssignVariableOp_12AssignVariableOp_122*
AssignVariableOp_13AssignVariableOp_132*
AssignVariableOp_14AssignVariableOp_142*
AssignVariableOp_15AssignVariableOp_152*
AssignVariableOp_16AssignVariableOp_162*
AssignVariableOp_17AssignVariableOp_172*
AssignVariableOp_18AssignVariableOp_182*
AssignVariableOp_19AssignVariableOp_192(
AssignVariableOp_2AssignVariableOp_22*
AssignVariableOp_20AssignVariableOp_202*
AssignVariableOp_21AssignVariableOp_212*
AssignVariableOp_22AssignVariableOp_222*
AssignVariableOp_23AssignVariableOp_232*
AssignVariableOp_24AssignVariableOp_242*
AssignVariableOp_25AssignVariableOp_252*
AssignVariableOp_26AssignVariableOp_262*
AssignVariableOp_27AssignVariableOp_272*
AssignVariableOp_28AssignVariableOp_282*
AssignVariableOp_29AssignVariableOp_292(
AssignVariableOp_3AssignVariableOp_32*
AssignVariableOp_30AssignVariableOp_302*
AssignVariableOp_31AssignVariableOp_312*
AssignVariableOp_32AssignVariableOp_322*
AssignVariableOp_33AssignVariableOp_332*
AssignVariableOp_34AssignVariableOp_342*
AssignVariableOp_35AssignVariableOp_352*
AssignVariableOp_36AssignVariableOp_362*
AssignVariableOp_37AssignVariableOp_372*
AssignVariableOp_38AssignVariableOp_382*
AssignVariableOp_39AssignVariableOp_392(
AssignVariableOp_4AssignVariableOp_42*
AssignVariableOp_40AssignVariableOp_402*
AssignVariableOp_41AssignVariableOp_412*
AssignVariableOp_42AssignVariableOp_422*
AssignVariableOp_43AssignVariableOp_432*
AssignVariableOp_44AssignVariableOp_442*
AssignVariableOp_45AssignVariableOp_452*
AssignVariableOp_46AssignVariableOp_462*
AssignVariableOp_47AssignVariableOp_472*
AssignVariableOp_48AssignVariableOp_482*
AssignVariableOp_49AssignVariableOp_492(
AssignVariableOp_5AssignVariableOp_52*
AssignVariableOp_50AssignVariableOp_502*
AssignVariableOp_51AssignVariableOp_512*
AssignVariableOp_52AssignVariableOp_522*
AssignVariableOp_53AssignVariableOp_532*
AssignVariableOp_54AssignVariableOp_542*
AssignVariableOp_55AssignVariableOp_552*
AssignVariableOp_56AssignVariableOp_562*
AssignVariableOp_57AssignVariableOp_572*
AssignVariableOp_58AssignVariableOp_582*
AssignVariableOp_59AssignVariableOp_592(
AssignVariableOp_6AssignVariableOp_62*
AssignVariableOp_60AssignVariableOp_602(
AssignVariableOp_7AssignVariableOp_72(
AssignVariableOp_8AssignVariableOp_82(
AssignVariableOp_9AssignVariableOp_9:C ?

_output_shapes
: 
%
_user_specified_namefile_prefix
�

�
'__inference_cvae_3_layer_call_fn_142361
input_1
unknown
	unknown_0
	unknown_1
	unknown_2
	unknown_3
	unknown_4
	unknown_5
	unknown_6
	unknown_7
	unknown_8
	unknown_9

unknown_10

unknown_11

unknown_12

unknown_13

unknown_14

unknown_15

unknown_16
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinput_1unknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6	unknown_7	unknown_8	unknown_9
unknown_10
unknown_11
unknown_12
unknown_13
unknown_14
unknown_15
unknown_16*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*4
_read_only_resource_inputs
	
*-
config_proto

CPU

GPU 2J 8� *K
fFRD
B__inference_cvae_3_layer_call_and_return_conditional_losses_1419642
StatefulPartitionedCall�
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*'
_output_shapes
:���������2

Identity"
identityIdentity:output:0*n
_input_shapes]
[:���������::::::::::::::::::22
StatefulPartitionedCallStatefulPartitionedCall:P L
'
_output_shapes
:���������
!
_user_specified_name	input_1
�
|
'__inference_z_mean_layer_call_fn_142990

inputs
unknown
	unknown_0
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *K
fFRD
B__inference_z_mean_layer_call_and_return_conditional_losses_1412472
StatefulPartitionedCall�
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*'
_output_shapes
:���������2

Identity"
identityIdentity:output:0*.
_input_shapes
:���������::22
StatefulPartitionedCallStatefulPartitionedCall:O K
'
_output_shapes
:���������
 
_user_specified_nameinputs
�
�
C__inference_decoder_layer_call_and_return_conditional_losses_141691

inputs
decoder_l2_141670
decoder_l2_141672
decoder_l3_141675
decoder_l3_141677
decoder_l4_141680
decoder_l4_141682
output_layer_141685
output_layer_141687
identity��"decoder_l2/StatefulPartitionedCall�"decoder_l3/StatefulPartitionedCall�"decoder_l4/StatefulPartitionedCall�$output_layer/StatefulPartitionedCall�
"decoder_l2/StatefulPartitionedCallStatefulPartitionedCallinputsdecoder_l2_141670decoder_l2_141672*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *O
fJRH
F__inference_decoder_l2_layer_call_and_return_conditional_losses_1414972$
"decoder_l2/StatefulPartitionedCall�
"decoder_l3/StatefulPartitionedCallStatefulPartitionedCall+decoder_l2/StatefulPartitionedCall:output:0decoder_l3_141675decoder_l3_141677*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *O
fJRH
F__inference_decoder_l3_layer_call_and_return_conditional_losses_1415242$
"decoder_l3/StatefulPartitionedCall�
"decoder_l4/StatefulPartitionedCallStatefulPartitionedCall+decoder_l3/StatefulPartitionedCall:output:0decoder_l4_141680decoder_l4_141682*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������
*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *O
fJRH
F__inference_decoder_l4_layer_call_and_return_conditional_losses_1415512$
"decoder_l4/StatefulPartitionedCall�
$output_layer/StatefulPartitionedCallStatefulPartitionedCall+decoder_l4/StatefulPartitionedCall:output:0output_layer_141685output_layer_141687*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *Q
fLRJ
H__inference_output_layer_layer_call_and_return_conditional_losses_1415782&
$output_layer/StatefulPartitionedCall�
IdentityIdentity-output_layer/StatefulPartitionedCall:output:0#^decoder_l2/StatefulPartitionedCall#^decoder_l3/StatefulPartitionedCall#^decoder_l4/StatefulPartitionedCall%^output_layer/StatefulPartitionedCall*
T0*'
_output_shapes
:���������2

Identity"
identityIdentity:output:0*F
_input_shapes5
3:���������::::::::2H
"decoder_l2/StatefulPartitionedCall"decoder_l2/StatefulPartitionedCall2H
"decoder_l3/StatefulPartitionedCall"decoder_l3/StatefulPartitionedCall2H
"decoder_l4/StatefulPartitionedCall"decoder_l4/StatefulPartitionedCall2L
$output_layer/StatefulPartitionedCall$output_layer/StatefulPartitionedCall:O K
'
_output_shapes
:���������
 
_user_specified_nameinputs"�L
saver_filename:0StatefulPartitionedCall_1:0StatefulPartitionedCall_28"
saved_model_main_op

NoOp*>
__saved_model_init_op%#
__saved_model_init_op

NoOp*�
serving_default�
;
input_10
serving_default_input_1:0���������<
output_10
StatefulPartitionedCall:0���������tensorflow/serving/predict:��
�
encoder
decoder
	optimizer
loss
	variables
trainable_variables
regularization_losses
	keras_api
	
signatures
�__call__
�_default_save_signature
+�&call_and_return_all_conditional_losses"�
_tf_keras_model�{"class_name": "CVAE", "name": "cvae_3", "trainable": true, "expects_training_arg": false, "dtype": "float32", "batch_input_shape": null, "must_restore_from_config": false, "config": {"layer was saved without config": true}, "is_graph_network": false, "keras_version": "2.4.0", "backend": "tensorflow", "model_config": {"class_name": "CVAE"}, "training_config": {"loss": null, "metrics": null, "weighted_metrics": null, "loss_weights": null, "optimizer_config": {"class_name": "Adam", "config": {"name": "Adam", "learning_rate": 0.009999999776482582, "decay": 0.0, "beta_1": 0.8999999761581421, "beta_2": 0.9990000128746033, "epsilon": 1e-07, "amsgrad": false}}}}
�7

layer-0
layer_with_weights-0
layer-1
layer_with_weights-1
layer-2
layer_with_weights-2
layer-3
layer_with_weights-3
layer-4
layer_with_weights-4
layer-5
layer-6
	variables
trainable_variables
regularization_losses
	keras_api
�__call__
+�&call_and_return_all_conditional_losses"�4
_tf_keras_network�4{"class_name": "Functional", "name": "encoder", "trainable": true, "expects_training_arg": true, "dtype": "float32", "batch_input_shape": null, "must_restore_from_config": false, "config": {"name": "encoder", "layers": [{"class_name": "InputLayer", "config": {"batch_input_shape": {"class_name": "__tuple__", "items": [null, 14]}, "dtype": "float32", "sparse": false, "ragged": false, "name": "encoder_input"}, "name": "encoder_input", "inbound_nodes": []}, {"class_name": "Dense", "config": {"name": "encoder_l2", "trainable": true, "dtype": "float32", "units": 10, "activation": "relu", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "name": "encoder_l2", "inbound_nodes": [[["encoder_input", 0, 0, {}]]]}, {"class_name": "Dense", "config": {"name": "encoder_l3", "trainable": true, "dtype": "float32", "units": 8, "activation": "relu", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "name": "encoder_l3", "inbound_nodes": [[["encoder_l2", 0, 0, {}]]]}, {"class_name": "Dense", "config": {"name": "encoder_l4", "trainable": true, "dtype": "float32", "units": 6, "activation": "relu", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "name": "encoder_l4", "inbound_nodes": [[["encoder_l3", 0, 0, {}]]]}, {"class_name": "Dense", "config": {"name": "z_mean", "trainable": true, "dtype": "float32", "units": 2, "activation": "linear", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "name": "z_mean", "inbound_nodes": [[["encoder_l4", 0, 0, {}]]]}, {"class_name": "Dense", "config": {"name": "z_log_var", "trainable": true, "dtype": "float32", "units": 2, "activation": "linear", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "name": "z_log_var", "inbound_nodes": [[["encoder_l4", 0, 0, {}]]]}, {"class_name": "Sampling", "config": {"name": "sampling_3", "trainable": true, "dtype": "float32"}, "name": "sampling_3", "inbound_nodes": [[["z_mean", 0, 0, {}], ["z_log_var", 0, 0, {}]]]}], "input_layers": [["encoder_input", 0, 0]], "output_layers": [["z_mean", 0, 0], ["z_log_var", 0, 0], ["sampling_3", 0, 0]]}, "input_spec": [{"class_name": "InputSpec", "config": {"dtype": null, "shape": {"class_name": "__tuple__", "items": [null, 14]}, "ndim": 2, "max_ndim": null, "min_ndim": null, "axes": {}}}], "build_input_shape": {"class_name": "TensorShape", "items": [null, 14]}, "is_graph_network": true, "keras_version": "2.4.0", "backend": "tensorflow", "model_config": {"class_name": "Functional", "config": {"name": "encoder", "layers": [{"class_name": "InputLayer", "config": {"batch_input_shape": {"class_name": "__tuple__", "items": [null, 14]}, "dtype": "float32", "sparse": false, "ragged": false, "name": "encoder_input"}, "name": "encoder_input", "inbound_nodes": []}, {"class_name": "Dense", "config": {"name": "encoder_l2", "trainable": true, "dtype": "float32", "units": 10, "activation": "relu", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "name": "encoder_l2", "inbound_nodes": [[["encoder_input", 0, 0, {}]]]}, {"class_name": "Dense", "config": {"name": "encoder_l3", "trainable": true, "dtype": "float32", "units": 8, "activation": "relu", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "name": "encoder_l3", "inbound_nodes": [[["encoder_l2", 0, 0, {}]]]}, {"class_name": "Dense", "config": {"name": "encoder_l4", "trainable": true, "dtype": "float32", "units": 6, "activation": "relu", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "name": "encoder_l4", "inbound_nodes": [[["encoder_l3", 0, 0, {}]]]}, {"class_name": "Dense", "config": {"name": "z_mean", "trainable": true, "dtype": "float32", "units": 2, "activation": "linear", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "name": "z_mean", "inbound_nodes": [[["encoder_l4", 0, 0, {}]]]}, {"class_name": "Dense", "config": {"name": "z_log_var", "trainable": true, "dtype": "float32", "units": 2, "activation": "linear", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "name": "z_log_var", "inbound_nodes": [[["encoder_l4", 0, 0, {}]]]}, {"class_name": "Sampling", "config": {"name": "sampling_3", "trainable": true, "dtype": "float32"}, "name": "sampling_3", "inbound_nodes": [[["z_mean", 0, 0, {}], ["z_log_var", 0, 0, {}]]]}], "input_layers": [["encoder_input", 0, 0]], "output_layers": [["z_mean", 0, 0], ["z_log_var", 0, 0], ["sampling_3", 0, 0]]}}}
�+
layer-0
layer_with_weights-0
layer-1
layer_with_weights-1
layer-2
layer_with_weights-2
layer-3
layer_with_weights-3
layer-4
	variables
trainable_variables
regularization_losses
	keras_api
�__call__
+�&call_and_return_all_conditional_losses"�)
_tf_keras_network�){"class_name": "Functional", "name": "decoder", "trainable": true, "expects_training_arg": true, "dtype": "float32", "batch_input_shape": null, "must_restore_from_config": false, "config": {"name": "decoder", "layers": [{"class_name": "InputLayer", "config": {"batch_input_shape": {"class_name": "__tuple__", "items": [null, 4]}, "dtype": "float32", "sparse": false, "ragged": false, "name": "latent_variable"}, "name": "latent_variable", "inbound_nodes": []}, {"class_name": "Dense", "config": {"name": "decoder_l2", "trainable": true, "dtype": "float32", "units": 6, "activation": "relu", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "name": "decoder_l2", "inbound_nodes": [[["latent_variable", 0, 0, {}]]]}, {"class_name": "Dense", "config": {"name": "decoder_l3", "trainable": true, "dtype": "float32", "units": 8, "activation": "relu", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "name": "decoder_l3", "inbound_nodes": [[["decoder_l2", 0, 0, {}]]]}, {"class_name": "Dense", "config": {"name": "decoder_l4", "trainable": true, "dtype": "float32", "units": 10, "activation": "relu", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "name": "decoder_l4", "inbound_nodes": [[["decoder_l3", 0, 0, {}]]]}, {"class_name": "Dense", "config": {"name": "output_layer", "trainable": true, "dtype": "float32", "units": 14, "activation": "relu", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "name": "output_layer", "inbound_nodes": [[["decoder_l4", 0, 0, {}]]]}], "input_layers": [["latent_variable", 0, 0]], "output_layers": [["output_layer", 0, 0]]}, "input_spec": [{"class_name": "InputSpec", "config": {"dtype": null, "shape": {"class_name": "__tuple__", "items": [null, 4]}, "ndim": 2, "max_ndim": null, "min_ndim": null, "axes": {}}}], "build_input_shape": {"class_name": "TensorShape", "items": [null, 4]}, "is_graph_network": true, "keras_version": "2.4.0", "backend": "tensorflow", "model_config": {"class_name": "Functional", "config": {"name": "decoder", "layers": [{"class_name": "InputLayer", "config": {"batch_input_shape": {"class_name": "__tuple__", "items": [null, 4]}, "dtype": "float32", "sparse": false, "ragged": false, "name": "latent_variable"}, "name": "latent_variable", "inbound_nodes": []}, {"class_name": "Dense", "config": {"name": "decoder_l2", "trainable": true, "dtype": "float32", "units": 6, "activation": "relu", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "name": "decoder_l2", "inbound_nodes": [[["latent_variable", 0, 0, {}]]]}, {"class_name": "Dense", "config": {"name": "decoder_l3", "trainable": true, "dtype": "float32", "units": 8, "activation": "relu", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "name": "decoder_l3", "inbound_nodes": [[["decoder_l2", 0, 0, {}]]]}, {"class_name": "Dense", "config": {"name": "decoder_l4", "trainable": true, "dtype": "float32", "units": 10, "activation": "relu", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "name": "decoder_l4", "inbound_nodes": [[["decoder_l3", 0, 0, {}]]]}, {"class_name": "Dense", "config": {"name": "output_layer", "trainable": true, "dtype": "float32", "units": 14, "activation": "relu", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "name": "output_layer", "inbound_nodes": [[["decoder_l4", 0, 0, {}]]]}], "input_layers": [["latent_variable", 0, 0]], "output_layers": [["output_layer", 0, 0]]}}}
�
iter

beta_1

 beta_2
	!decay
"learning_rate#m�$m�%m�&m�'m�(m�)m�*m�+m�,m�-m�.m�/m�0m�1m�2m�3m�4m�#v�$v�%v�&v�'v�(v�)v�*v�+v�,v�-v�.v�/v�0v�1v�2v�3v�4v�"
	optimizer
 "
trackable_dict_wrapper
�
#0
$1
%2
&3
'4
(5
)6
*7
+8
,9
-10
.11
/12
013
114
215
316
417"
trackable_list_wrapper
�
#0
$1
%2
&3
'4
(5
)6
*7
+8
,9
-10
.11
/12
013
114
215
316
417"
trackable_list_wrapper
 "
trackable_list_wrapper
�
	variables
5metrics
6layer_regularization_losses
trainable_variables
7non_trainable_variables

8layers
9layer_metrics
regularization_losses
�__call__
�_default_save_signature
+�&call_and_return_all_conditional_losses
'�"call_and_return_conditional_losses"
_generic_user_object
-
�serving_default"
signature_map
�"�
_tf_keras_input_layer�{"class_name": "InputLayer", "name": "encoder_input", "dtype": "float32", "sparse": false, "ragged": false, "batch_input_shape": {"class_name": "__tuple__", "items": [null, 14]}, "config": {"batch_input_shape": {"class_name": "__tuple__", "items": [null, 14]}, "dtype": "float32", "sparse": false, "ragged": false, "name": "encoder_input"}}
�

#kernel
$bias
:	variables
;trainable_variables
<regularization_losses
=	keras_api
�__call__
+�&call_and_return_all_conditional_losses"�
_tf_keras_layer�{"class_name": "Dense", "name": "encoder_l2", "trainable": true, "expects_training_arg": false, "dtype": "float32", "batch_input_shape": null, "stateful": false, "must_restore_from_config": false, "config": {"name": "encoder_l2", "trainable": true, "dtype": "float32", "units": 10, "activation": "relu", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "input_spec": {"class_name": "InputSpec", "config": {"dtype": null, "shape": null, "ndim": null, "max_ndim": null, "min_ndim": 2, "axes": {"-1": 14}}}, "build_input_shape": {"class_name": "TensorShape", "items": [null, 14]}}
�

%kernel
&bias
>	variables
?trainable_variables
@regularization_losses
A	keras_api
�__call__
+�&call_and_return_all_conditional_losses"�
_tf_keras_layer�{"class_name": "Dense", "name": "encoder_l3", "trainable": true, "expects_training_arg": false, "dtype": "float32", "batch_input_shape": null, "stateful": false, "must_restore_from_config": false, "config": {"name": "encoder_l3", "trainable": true, "dtype": "float32", "units": 8, "activation": "relu", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "input_spec": {"class_name": "InputSpec", "config": {"dtype": null, "shape": null, "ndim": null, "max_ndim": null, "min_ndim": 2, "axes": {"-1": 10}}}, "build_input_shape": {"class_name": "TensorShape", "items": [null, 10]}}
�

'kernel
(bias
B	variables
Ctrainable_variables
Dregularization_losses
E	keras_api
�__call__
+�&call_and_return_all_conditional_losses"�
_tf_keras_layer�{"class_name": "Dense", "name": "encoder_l4", "trainable": true, "expects_training_arg": false, "dtype": "float32", "batch_input_shape": null, "stateful": false, "must_restore_from_config": false, "config": {"name": "encoder_l4", "trainable": true, "dtype": "float32", "units": 6, "activation": "relu", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "input_spec": {"class_name": "InputSpec", "config": {"dtype": null, "shape": null, "ndim": null, "max_ndim": null, "min_ndim": 2, "axes": {"-1": 8}}}, "build_input_shape": {"class_name": "TensorShape", "items": [null, 8]}}
�

)kernel
*bias
F	variables
Gtrainable_variables
Hregularization_losses
I	keras_api
�__call__
+�&call_and_return_all_conditional_losses"�
_tf_keras_layer�{"class_name": "Dense", "name": "z_mean", "trainable": true, "expects_training_arg": false, "dtype": "float32", "batch_input_shape": null, "stateful": false, "must_restore_from_config": false, "config": {"name": "z_mean", "trainable": true, "dtype": "float32", "units": 2, "activation": "linear", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "input_spec": {"class_name": "InputSpec", "config": {"dtype": null, "shape": null, "ndim": null, "max_ndim": null, "min_ndim": 2, "axes": {"-1": 6}}}, "build_input_shape": {"class_name": "TensorShape", "items": [null, 6]}}
�

+kernel
,bias
J	variables
Ktrainable_variables
Lregularization_losses
M	keras_api
�__call__
+�&call_and_return_all_conditional_losses"�
_tf_keras_layer�{"class_name": "Dense", "name": "z_log_var", "trainable": true, "expects_training_arg": false, "dtype": "float32", "batch_input_shape": null, "stateful": false, "must_restore_from_config": false, "config": {"name": "z_log_var", "trainable": true, "dtype": "float32", "units": 2, "activation": "linear", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "input_spec": {"class_name": "InputSpec", "config": {"dtype": null, "shape": null, "ndim": null, "max_ndim": null, "min_ndim": 2, "axes": {"-1": 6}}}, "build_input_shape": {"class_name": "TensorShape", "items": [null, 6]}}
�
N	variables
Otrainable_variables
Pregularization_losses
Q	keras_api
�__call__
+�&call_and_return_all_conditional_losses"�
_tf_keras_layer�{"class_name": "Sampling", "name": "sampling_3", "trainable": true, "expects_training_arg": false, "dtype": "float32", "batch_input_shape": null, "stateful": false, "must_restore_from_config": false, "config": {"name": "sampling_3", "trainable": true, "dtype": "float32"}}
f
#0
$1
%2
&3
'4
(5
)6
*7
+8
,9"
trackable_list_wrapper
f
#0
$1
%2
&3
'4
(5
)6
*7
+8
,9"
trackable_list_wrapper
 "
trackable_list_wrapper
�
	variables
Rmetrics
Slayer_regularization_losses
trainable_variables
Tnon_trainable_variables

Ulayers
Vlayer_metrics
regularization_losses
�__call__
+�&call_and_return_all_conditional_losses
'�"call_and_return_conditional_losses"
_generic_user_object
�"�
_tf_keras_input_layer�{"class_name": "InputLayer", "name": "latent_variable", "dtype": "float32", "sparse": false, "ragged": false, "batch_input_shape": {"class_name": "__tuple__", "items": [null, 4]}, "config": {"batch_input_shape": {"class_name": "__tuple__", "items": [null, 4]}, "dtype": "float32", "sparse": false, "ragged": false, "name": "latent_variable"}}
�

-kernel
.bias
W	variables
Xtrainable_variables
Yregularization_losses
Z	keras_api
�__call__
+�&call_and_return_all_conditional_losses"�
_tf_keras_layer�{"class_name": "Dense", "name": "decoder_l2", "trainable": true, "expects_training_arg": false, "dtype": "float32", "batch_input_shape": null, "stateful": false, "must_restore_from_config": false, "config": {"name": "decoder_l2", "trainable": true, "dtype": "float32", "units": 6, "activation": "relu", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "input_spec": {"class_name": "InputSpec", "config": {"dtype": null, "shape": null, "ndim": null, "max_ndim": null, "min_ndim": 2, "axes": {"-1": 4}}}, "build_input_shape": {"class_name": "TensorShape", "items": [null, 4]}}
�

/kernel
0bias
[	variables
\trainable_variables
]regularization_losses
^	keras_api
�__call__
+�&call_and_return_all_conditional_losses"�
_tf_keras_layer�{"class_name": "Dense", "name": "decoder_l3", "trainable": true, "expects_training_arg": false, "dtype": "float32", "batch_input_shape": null, "stateful": false, "must_restore_from_config": false, "config": {"name": "decoder_l3", "trainable": true, "dtype": "float32", "units": 8, "activation": "relu", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "input_spec": {"class_name": "InputSpec", "config": {"dtype": null, "shape": null, "ndim": null, "max_ndim": null, "min_ndim": 2, "axes": {"-1": 6}}}, "build_input_shape": {"class_name": "TensorShape", "items": [null, 6]}}
�

1kernel
2bias
_	variables
`trainable_variables
aregularization_losses
b	keras_api
�__call__
+�&call_and_return_all_conditional_losses"�
_tf_keras_layer�{"class_name": "Dense", "name": "decoder_l4", "trainable": true, "expects_training_arg": false, "dtype": "float32", "batch_input_shape": null, "stateful": false, "must_restore_from_config": false, "config": {"name": "decoder_l4", "trainable": true, "dtype": "float32", "units": 10, "activation": "relu", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "input_spec": {"class_name": "InputSpec", "config": {"dtype": null, "shape": null, "ndim": null, "max_ndim": null, "min_ndim": 2, "axes": {"-1": 8}}}, "build_input_shape": {"class_name": "TensorShape", "items": [null, 8]}}
�

3kernel
4bias
c	variables
dtrainable_variables
eregularization_losses
f	keras_api
�__call__
+�&call_and_return_all_conditional_losses"�
_tf_keras_layer�{"class_name": "Dense", "name": "output_layer", "trainable": true, "expects_training_arg": false, "dtype": "float32", "batch_input_shape": null, "stateful": false, "must_restore_from_config": false, "config": {"name": "output_layer", "trainable": true, "dtype": "float32", "units": 14, "activation": "relu", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "input_spec": {"class_name": "InputSpec", "config": {"dtype": null, "shape": null, "ndim": null, "max_ndim": null, "min_ndim": 2, "axes": {"-1": 10}}}, "build_input_shape": {"class_name": "TensorShape", "items": [null, 10]}}
X
-0
.1
/2
03
14
25
36
47"
trackable_list_wrapper
X
-0
.1
/2
03
14
25
36
47"
trackable_list_wrapper
 "
trackable_list_wrapper
�
	variables
gmetrics
hlayer_regularization_losses
trainable_variables
inon_trainable_variables

jlayers
klayer_metrics
regularization_losses
�__call__
+�&call_and_return_all_conditional_losses
'�"call_and_return_conditional_losses"
_generic_user_object
:	 (2	Adam/iter
: (2Adam/beta_1
: (2Adam/beta_2
: (2
Adam/decay
: (2Adam/learning_rate
#:!
2encoder_l2/kernel
:
2encoder_l2/bias
#:!
2encoder_l3/kernel
:2encoder_l3/bias
#:!2encoder_l4/kernel
:2encoder_l4/bias
:2z_mean/kernel
:2z_mean/bias
": 2z_log_var/kernel
:2z_log_var/bias
#:!2decoder_l2/kernel
:2decoder_l2/bias
#:!2decoder_l3/kernel
:2decoder_l3/bias
#:!
2decoder_l4/kernel
:
2decoder_l4/bias
%:#
2output_layer/kernel
:2output_layer/bias
'
l0"
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
.
0
1"
trackable_list_wrapper
 "
trackable_dict_wrapper
.
#0
$1"
trackable_list_wrapper
.
#0
$1"
trackable_list_wrapper
 "
trackable_list_wrapper
�
:	variables
mmetrics
nlayer_regularization_losses
;trainable_variables
onon_trainable_variables

players
qlayer_metrics
<regularization_losses
�__call__
+�&call_and_return_all_conditional_losses
'�"call_and_return_conditional_losses"
_generic_user_object
.
%0
&1"
trackable_list_wrapper
.
%0
&1"
trackable_list_wrapper
 "
trackable_list_wrapper
�
>	variables
rmetrics
slayer_regularization_losses
?trainable_variables
tnon_trainable_variables

ulayers
vlayer_metrics
@regularization_losses
�__call__
+�&call_and_return_all_conditional_losses
'�"call_and_return_conditional_losses"
_generic_user_object
.
'0
(1"
trackable_list_wrapper
.
'0
(1"
trackable_list_wrapper
 "
trackable_list_wrapper
�
B	variables
wmetrics
xlayer_regularization_losses
Ctrainable_variables
ynon_trainable_variables

zlayers
{layer_metrics
Dregularization_losses
�__call__
+�&call_and_return_all_conditional_losses
'�"call_and_return_conditional_losses"
_generic_user_object
.
)0
*1"
trackable_list_wrapper
.
)0
*1"
trackable_list_wrapper
 "
trackable_list_wrapper
�
F	variables
|metrics
}layer_regularization_losses
Gtrainable_variables
~non_trainable_variables

layers
�layer_metrics
Hregularization_losses
�__call__
+�&call_and_return_all_conditional_losses
'�"call_and_return_conditional_losses"
_generic_user_object
.
+0
,1"
trackable_list_wrapper
.
+0
,1"
trackable_list_wrapper
 "
trackable_list_wrapper
�
J	variables
�metrics
 �layer_regularization_losses
Ktrainable_variables
�non_trainable_variables
�layers
�layer_metrics
Lregularization_losses
�__call__
+�&call_and_return_all_conditional_losses
'�"call_and_return_conditional_losses"
_generic_user_object
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
�
N	variables
�metrics
 �layer_regularization_losses
Otrainable_variables
�non_trainable_variables
�layers
�layer_metrics
Pregularization_losses
�__call__
+�&call_and_return_all_conditional_losses
'�"call_and_return_conditional_losses"
_generic_user_object
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
Q

0
1
2
3
4
5
6"
trackable_list_wrapper
 "
trackable_dict_wrapper
.
-0
.1"
trackable_list_wrapper
.
-0
.1"
trackable_list_wrapper
 "
trackable_list_wrapper
�
W	variables
�metrics
 �layer_regularization_losses
Xtrainable_variables
�non_trainable_variables
�layers
�layer_metrics
Yregularization_losses
�__call__
+�&call_and_return_all_conditional_losses
'�"call_and_return_conditional_losses"
_generic_user_object
.
/0
01"
trackable_list_wrapper
.
/0
01"
trackable_list_wrapper
 "
trackable_list_wrapper
�
[	variables
�metrics
 �layer_regularization_losses
\trainable_variables
�non_trainable_variables
�layers
�layer_metrics
]regularization_losses
�__call__
+�&call_and_return_all_conditional_losses
'�"call_and_return_conditional_losses"
_generic_user_object
.
10
21"
trackable_list_wrapper
.
10
21"
trackable_list_wrapper
 "
trackable_list_wrapper
�
_	variables
�metrics
 �layer_regularization_losses
`trainable_variables
�non_trainable_variables
�layers
�layer_metrics
aregularization_losses
�__call__
+�&call_and_return_all_conditional_losses
'�"call_and_return_conditional_losses"
_generic_user_object
.
30
41"
trackable_list_wrapper
.
30
41"
trackable_list_wrapper
 "
trackable_list_wrapper
�
c	variables
�metrics
 �layer_regularization_losses
dtrainable_variables
�non_trainable_variables
�layers
�layer_metrics
eregularization_losses
�__call__
+�&call_and_return_all_conditional_losses
'�"call_and_return_conditional_losses"
_generic_user_object
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
C
0
1
2
3
4"
trackable_list_wrapper
 "
trackable_dict_wrapper
�

�total

�count
�	variables
�	keras_api"�
_tf_keras_metricj{"class_name": "Mean", "name": "loss", "dtype": "float32", "config": {"name": "loss", "dtype": "float32"}}
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
:  (2total
:  (2count
0
�0
�1"
trackable_list_wrapper
.
�	variables"
_generic_user_object
(:&
2Adam/encoder_l2/kernel/m
": 
2Adam/encoder_l2/bias/m
(:&
2Adam/encoder_l3/kernel/m
": 2Adam/encoder_l3/bias/m
(:&2Adam/encoder_l4/kernel/m
": 2Adam/encoder_l4/bias/m
$:"2Adam/z_mean/kernel/m
:2Adam/z_mean/bias/m
':%2Adam/z_log_var/kernel/m
!:2Adam/z_log_var/bias/m
(:&2Adam/decoder_l2/kernel/m
": 2Adam/decoder_l2/bias/m
(:&2Adam/decoder_l3/kernel/m
": 2Adam/decoder_l3/bias/m
(:&
2Adam/decoder_l4/kernel/m
": 
2Adam/decoder_l4/bias/m
*:(
2Adam/output_layer/kernel/m
$:"2Adam/output_layer/bias/m
(:&
2Adam/encoder_l2/kernel/v
": 
2Adam/encoder_l2/bias/v
(:&
2Adam/encoder_l3/kernel/v
": 2Adam/encoder_l3/bias/v
(:&2Adam/encoder_l4/kernel/v
": 2Adam/encoder_l4/bias/v
$:"2Adam/z_mean/kernel/v
:2Adam/z_mean/bias/v
':%2Adam/z_log_var/kernel/v
!:2Adam/z_log_var/bias/v
(:&2Adam/decoder_l2/kernel/v
": 2Adam/decoder_l2/bias/v
(:&2Adam/decoder_l3/kernel/v
": 2Adam/decoder_l3/bias/v
(:&
2Adam/decoder_l4/kernel/v
": 
2Adam/decoder_l4/bias/v
*:(
2Adam/output_layer/kernel/v
$:"2Adam/output_layer/bias/v
�2�
'__inference_cvae_3_layer_call_fn_142361
'__inference_cvae_3_layer_call_fn_142320
'__inference_cvae_3_layer_call_fn_142586
'__inference_cvae_3_layer_call_fn_142627�
���
FullArgSpec)
args!�
jself
jinputs

jtraining
varargs
 
varkw
 
defaults�
p 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�2�
!__inference__wrapped_model_141152�
���
FullArgSpec
args� 
varargsjargs
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *&�#
!�
input_1���������
�2�
B__inference_cvae_3_layer_call_and_return_conditional_losses_142187
B__inference_cvae_3_layer_call_and_return_conditional_losses_142279
B__inference_cvae_3_layer_call_and_return_conditional_losses_142453
B__inference_cvae_3_layer_call_and_return_conditional_losses_142545�
���
FullArgSpec)
args!�
jself
jinputs

jtraining
varargs
 
varkw
 
defaults�
p 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�2�
(__inference_encoder_layer_call_fn_141421
(__inference_encoder_layer_call_fn_142805
(__inference_encoder_layer_call_fn_141482
(__inference_encoder_layer_call_fn_142776�
���
FullArgSpec1
args)�&
jself
jinputs

jtraining
jmask
varargs
 
varkw
 
defaults�
p 

 

kwonlyargs� 
kwonlydefaults� 
annotations� *
 
�2�
C__inference_encoder_layer_call_and_return_conditional_losses_142747
C__inference_encoder_layer_call_and_return_conditional_losses_141359
C__inference_encoder_layer_call_and_return_conditional_losses_141327
C__inference_encoder_layer_call_and_return_conditional_losses_142687�
���
FullArgSpec1
args)�&
jself
jinputs

jtraining
jmask
varargs
 
varkw
 
defaults�
p 

 

kwonlyargs� 
kwonlydefaults� 
annotations� *
 
�2�
(__inference_decoder_layer_call_fn_141665
(__inference_decoder_layer_call_fn_141710
(__inference_decoder_layer_call_fn_142890
(__inference_decoder_layer_call_fn_142911�
���
FullArgSpec1
args)�&
jself
jinputs

jtraining
jmask
varargs
 
varkw
 
defaults�
p 

 

kwonlyargs� 
kwonlydefaults� 
annotations� *
 
�2�
C__inference_decoder_layer_call_and_return_conditional_losses_142869
C__inference_decoder_layer_call_and_return_conditional_losses_142837
C__inference_decoder_layer_call_and_return_conditional_losses_141595
C__inference_decoder_layer_call_and_return_conditional_losses_141619�
���
FullArgSpec1
args)�&
jself
jinputs

jtraining
jmask
varargs
 
varkw
 
defaults�
p 

 

kwonlyargs� 
kwonlydefaults� 
annotations� *
 
�B�
$__inference_signature_wrapper_142095input_1"�
���
FullArgSpec
args� 
varargs
 
varkwjkwargs
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�2�
+__inference_encoder_l2_layer_call_fn_142931�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�2�
F__inference_encoder_l2_layer_call_and_return_conditional_losses_142922�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�2�
+__inference_encoder_l3_layer_call_fn_142951�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�2�
F__inference_encoder_l3_layer_call_and_return_conditional_losses_142942�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�2�
+__inference_encoder_l4_layer_call_fn_142971�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�2�
F__inference_encoder_l4_layer_call_and_return_conditional_losses_142962�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�2�
'__inference_z_mean_layer_call_fn_142990�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�2�
B__inference_z_mean_layer_call_and_return_conditional_losses_142981�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�2�
*__inference_z_log_var_layer_call_fn_143009�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�2�
E__inference_z_log_var_layer_call_and_return_conditional_losses_143000�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�2�
+__inference_sampling_3_layer_call_fn_143041�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�2�
F__inference_sampling_3_layer_call_and_return_conditional_losses_143035�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�2�
+__inference_decoder_l2_layer_call_fn_143061�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�2�
F__inference_decoder_l2_layer_call_and_return_conditional_losses_143052�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�2�
+__inference_decoder_l3_layer_call_fn_143081�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�2�
F__inference_decoder_l3_layer_call_and_return_conditional_losses_143072�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�2�
+__inference_decoder_l4_layer_call_fn_143101�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�2�
F__inference_decoder_l4_layer_call_and_return_conditional_losses_143092�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�2�
-__inference_output_layer_layer_call_fn_143121�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�2�
H__inference_output_layer_layer_call_and_return_conditional_losses_143112�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 �
!__inference__wrapped_model_141152{#$%&'()*+,-./012340�-
&�#
!�
input_1���������
� "3�0
.
output_1"�
output_1����������
B__inference_cvae_3_layer_call_and_return_conditional_losses_142187q#$%&'()*+,-./012344�1
*�'
!�
input_1���������
p
� "%�"
�
0���������
� �
B__inference_cvae_3_layer_call_and_return_conditional_losses_142279q#$%&'()*+,-./012344�1
*�'
!�
input_1���������
p 
� "%�"
�
0���������
� �
B__inference_cvae_3_layer_call_and_return_conditional_losses_142453p#$%&'()*+,-./012343�0
)�&
 �
inputs���������
p
� "%�"
�
0���������
� �
B__inference_cvae_3_layer_call_and_return_conditional_losses_142545p#$%&'()*+,-./012343�0
)�&
 �
inputs���������
p 
� "%�"
�
0���������
� �
'__inference_cvae_3_layer_call_fn_142320d#$%&'()*+,-./012344�1
*�'
!�
input_1���������
p
� "�����������
'__inference_cvae_3_layer_call_fn_142361d#$%&'()*+,-./012344�1
*�'
!�
input_1���������
p 
� "�����������
'__inference_cvae_3_layer_call_fn_142586c#$%&'()*+,-./012343�0
)�&
 �
inputs���������
p
� "�����������
'__inference_cvae_3_layer_call_fn_142627c#$%&'()*+,-./012343�0
)�&
 �
inputs���������
p 
� "�����������
F__inference_decoder_l2_layer_call_and_return_conditional_losses_143052\-./�,
%�"
 �
inputs���������
� "%�"
�
0���������
� ~
+__inference_decoder_l2_layer_call_fn_143061O-./�,
%�"
 �
inputs���������
� "�����������
F__inference_decoder_l3_layer_call_and_return_conditional_losses_143072\/0/�,
%�"
 �
inputs���������
� "%�"
�
0���������
� ~
+__inference_decoder_l3_layer_call_fn_143081O/0/�,
%�"
 �
inputs���������
� "�����������
F__inference_decoder_l4_layer_call_and_return_conditional_losses_143092\12/�,
%�"
 �
inputs���������
� "%�"
�
0���������

� ~
+__inference_decoder_l4_layer_call_fn_143101O12/�,
%�"
 �
inputs���������
� "����������
�
C__inference_decoder_layer_call_and_return_conditional_losses_141595s-./01234@�=
6�3
)�&
latent_variable���������
p

 
� "%�"
�
0���������
� �
C__inference_decoder_layer_call_and_return_conditional_losses_141619s-./01234@�=
6�3
)�&
latent_variable���������
p 

 
� "%�"
�
0���������
� �
C__inference_decoder_layer_call_and_return_conditional_losses_142837j-./012347�4
-�*
 �
inputs���������
p

 
� "%�"
�
0���������
� �
C__inference_decoder_layer_call_and_return_conditional_losses_142869j-./012347�4
-�*
 �
inputs���������
p 

 
� "%�"
�
0���������
� �
(__inference_decoder_layer_call_fn_141665f-./01234@�=
6�3
)�&
latent_variable���������
p

 
� "�����������
(__inference_decoder_layer_call_fn_141710f-./01234@�=
6�3
)�&
latent_variable���������
p 

 
� "�����������
(__inference_decoder_layer_call_fn_142890]-./012347�4
-�*
 �
inputs���������
p

 
� "�����������
(__inference_decoder_layer_call_fn_142911]-./012347�4
-�*
 �
inputs���������
p 

 
� "�����������
F__inference_encoder_l2_layer_call_and_return_conditional_losses_142922\#$/�,
%�"
 �
inputs���������
� "%�"
�
0���������

� ~
+__inference_encoder_l2_layer_call_fn_142931O#$/�,
%�"
 �
inputs���������
� "����������
�
F__inference_encoder_l3_layer_call_and_return_conditional_losses_142942\%&/�,
%�"
 �
inputs���������

� "%�"
�
0���������
� ~
+__inference_encoder_l3_layer_call_fn_142951O%&/�,
%�"
 �
inputs���������

� "�����������
F__inference_encoder_l4_layer_call_and_return_conditional_losses_142962\'(/�,
%�"
 �
inputs���������
� "%�"
�
0���������
� ~
+__inference_encoder_l4_layer_call_fn_142971O'(/�,
%�"
 �
inputs���������
� "�����������
C__inference_encoder_layer_call_and_return_conditional_losses_141327�
#$%&'()*+,>�;
4�1
'�$
encoder_input���������
p

 
� "j�g
`�]
�
0/0���������
�
0/1���������
�
0/2���������
� �
C__inference_encoder_layer_call_and_return_conditional_losses_141359�
#$%&'()*+,>�;
4�1
'�$
encoder_input���������
p 

 
� "j�g
`�]
�
0/0���������
�
0/1���������
�
0/2���������
� �
C__inference_encoder_layer_call_and_return_conditional_losses_142687�
#$%&'()*+,7�4
-�*
 �
inputs���������
p

 
� "j�g
`�]
�
0/0���������
�
0/1���������
�
0/2���������
� �
C__inference_encoder_layer_call_and_return_conditional_losses_142747�
#$%&'()*+,7�4
-�*
 �
inputs���������
p 

 
� "j�g
`�]
�
0/0���������
�
0/1���������
�
0/2���������
� �
(__inference_encoder_layer_call_fn_141421�
#$%&'()*+,>�;
4�1
'�$
encoder_input���������
p

 
� "Z�W
�
0���������
�
1���������
�
2����������
(__inference_encoder_layer_call_fn_141482�
#$%&'()*+,>�;
4�1
'�$
encoder_input���������
p 

 
� "Z�W
�
0���������
�
1���������
�
2����������
(__inference_encoder_layer_call_fn_142776�
#$%&'()*+,7�4
-�*
 �
inputs���������
p

 
� "Z�W
�
0���������
�
1���������
�
2����������
(__inference_encoder_layer_call_fn_142805�
#$%&'()*+,7�4
-�*
 �
inputs���������
p 

 
� "Z�W
�
0���������
�
1���������
�
2����������
H__inference_output_layer_layer_call_and_return_conditional_losses_143112\34/�,
%�"
 �
inputs���������

� "%�"
�
0���������
� �
-__inference_output_layer_layer_call_fn_143121O34/�,
%�"
 �
inputs���������

� "�����������
F__inference_sampling_3_layer_call_and_return_conditional_losses_143035�Z�W
P�M
K�H
"�
inputs/0���������
"�
inputs/1���������
� "%�"
�
0���������
� �
+__inference_sampling_3_layer_call_fn_143041vZ�W
P�M
K�H
"�
inputs/0���������
"�
inputs/1���������
� "�����������
$__inference_signature_wrapper_142095�#$%&'()*+,-./01234;�8
� 
1�.
,
input_1!�
input_1���������"3�0
.
output_1"�
output_1����������
E__inference_z_log_var_layer_call_and_return_conditional_losses_143000\+,/�,
%�"
 �
inputs���������
� "%�"
�
0���������
� }
*__inference_z_log_var_layer_call_fn_143009O+,/�,
%�"
 �
inputs���������
� "�����������
B__inference_z_mean_layer_call_and_return_conditional_losses_142981\)*/�,
%�"
 �
inputs���������
� "%�"
�
0���������
� z
'__inference_z_mean_layer_call_fn_142990O)*/�,
%�"
 �
inputs���������
� "����������