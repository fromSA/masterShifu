лФ
фЕ
B
AssignVariableOp
resource
value"dtype"
dtypetype
~
BiasAdd

value"T	
bias"T
output"T" 
Ttype:
2	"-
data_formatstringNHWC:
NHWCNCHW
8
Const
output"dtype"
valuetensor"
dtypetype
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
delete_old_dirsbool(
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
@
ReadVariableOp
resource
value"dtype"
dtypetype
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
list(type)(0
l
SaveV2

prefix
tensor_names
shape_and_slices
tensors2dtypes"
dtypes
list(type)(0
?
Select
	condition

t"T
e"T
output"T"	
Ttype
H
ShardedFilename
basename	
shard

num_shards
filename
0
Sigmoid
x"T
y"T"
Ttype:

2
О
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
executor_typestring 
@
StaticRegexFullMatch	
input

output
"
patternstring
N

StringJoin
inputs*N

output"
Nint(0"
	separatorstring 

VarHandleOp
resource"
	containerstring "
shared_namestring "
dtypetype"
shapeshape"#
allowed_deviceslist(string)
 "serve*2.4.02v2.4.0-rc4-71-g582c8d236cb8Лг
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

generator_l1/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:
*$
shared_namegenerator_l1/kernel
{
'generator_l1/kernel/Read/ReadVariableOpReadVariableOpgenerator_l1/kernel*
_output_shapes

:
*
dtype0
z
generator_l1/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:
*"
shared_namegenerator_l1/bias
s
%generator_l1/bias/Read/ReadVariableOpReadVariableOpgenerator_l1/bias*
_output_shapes
:
*
dtype0

generator_l2/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:
*$
shared_namegenerator_l2/kernel
{
'generator_l2/kernel/Read/ReadVariableOpReadVariableOpgenerator_l2/kernel*
_output_shapes

:
*
dtype0
z
generator_l2/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:*"
shared_namegenerator_l2/bias
s
%generator_l2/bias/Read/ReadVariableOpReadVariableOpgenerator_l2/bias*
_output_shapes
:*
dtype0

generator_l3/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:*$
shared_namegenerator_l3/kernel
{
'generator_l3/kernel/Read/ReadVariableOpReadVariableOpgenerator_l3/kernel*
_output_shapes

:*
dtype0
z
generator_l3/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:*"
shared_namegenerator_l3/bias
s
%generator_l3/bias/Read/ReadVariableOpReadVariableOpgenerator_l3/bias*
_output_shapes
:*
dtype0

generator_output/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:*(
shared_namegenerator_output/kernel

+generator_output/kernel/Read/ReadVariableOpReadVariableOpgenerator_output/kernel*
_output_shapes

:*
dtype0

generator_output/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:*&
shared_namegenerator_output/bias
{
)generator_output/bias/Read/ReadVariableOpReadVariableOpgenerator_output/bias*
_output_shapes
:*
dtype0

discriminator_l1/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:
*(
shared_namediscriminator_l1/kernel

+discriminator_l1/kernel/Read/ReadVariableOpReadVariableOpdiscriminator_l1/kernel*
_output_shapes

:
*
dtype0

discriminator_l1/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:
*&
shared_namediscriminator_l1/bias
{
)discriminator_l1/bias/Read/ReadVariableOpReadVariableOpdiscriminator_l1/bias*
_output_shapes
:
*
dtype0

discriminator_l2/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:
*(
shared_namediscriminator_l2/kernel

+discriminator_l2/kernel/Read/ReadVariableOpReadVariableOpdiscriminator_l2/kernel*
_output_shapes

:
*
dtype0

discriminator_l2/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:*&
shared_namediscriminator_l2/bias
{
)discriminator_l2/bias/Read/ReadVariableOpReadVariableOpdiscriminator_l2/bias*
_output_shapes
:*
dtype0

discriminator_l3/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:*(
shared_namediscriminator_l3/kernel

+discriminator_l3/kernel/Read/ReadVariableOpReadVariableOpdiscriminator_l3/kernel*
_output_shapes

:*
dtype0

discriminator_l3/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:*&
shared_namediscriminator_l3/bias
{
)discriminator_l3/bias/Read/ReadVariableOpReadVariableOpdiscriminator_l3/bias*
_output_shapes
:*
dtype0

discriminator_output/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:*,
shared_namediscriminator_output/kernel

/discriminator_output/kernel/Read/ReadVariableOpReadVariableOpdiscriminator_output/kernel*
_output_shapes

:*
dtype0

discriminator_output/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:**
shared_namediscriminator_output/bias

-discriminator_output/bias/Read/ReadVariableOpReadVariableOpdiscriminator_output/bias*
_output_shapes
:*
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

Adam/generator_l1/kernel/mVarHandleOp*
_output_shapes
: *
dtype0*
shape
:
*+
shared_nameAdam/generator_l1/kernel/m

.Adam/generator_l1/kernel/m/Read/ReadVariableOpReadVariableOpAdam/generator_l1/kernel/m*
_output_shapes

:
*
dtype0

Adam/generator_l1/bias/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:
*)
shared_nameAdam/generator_l1/bias/m

,Adam/generator_l1/bias/m/Read/ReadVariableOpReadVariableOpAdam/generator_l1/bias/m*
_output_shapes
:
*
dtype0

Adam/generator_l2/kernel/mVarHandleOp*
_output_shapes
: *
dtype0*
shape
:
*+
shared_nameAdam/generator_l2/kernel/m

.Adam/generator_l2/kernel/m/Read/ReadVariableOpReadVariableOpAdam/generator_l2/kernel/m*
_output_shapes

:
*
dtype0

Adam/generator_l2/bias/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:*)
shared_nameAdam/generator_l2/bias/m

,Adam/generator_l2/bias/m/Read/ReadVariableOpReadVariableOpAdam/generator_l2/bias/m*
_output_shapes
:*
dtype0

Adam/generator_l3/kernel/mVarHandleOp*
_output_shapes
: *
dtype0*
shape
:*+
shared_nameAdam/generator_l3/kernel/m

.Adam/generator_l3/kernel/m/Read/ReadVariableOpReadVariableOpAdam/generator_l3/kernel/m*
_output_shapes

:*
dtype0

Adam/generator_l3/bias/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:*)
shared_nameAdam/generator_l3/bias/m

,Adam/generator_l3/bias/m/Read/ReadVariableOpReadVariableOpAdam/generator_l3/bias/m*
_output_shapes
:*
dtype0

Adam/generator_output/kernel/mVarHandleOp*
_output_shapes
: *
dtype0*
shape
:*/
shared_name Adam/generator_output/kernel/m

2Adam/generator_output/kernel/m/Read/ReadVariableOpReadVariableOpAdam/generator_output/kernel/m*
_output_shapes

:*
dtype0

Adam/generator_output/bias/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:*-
shared_nameAdam/generator_output/bias/m

0Adam/generator_output/bias/m/Read/ReadVariableOpReadVariableOpAdam/generator_output/bias/m*
_output_shapes
:*
dtype0

 Adam_1/discriminator_l1/kernel/mVarHandleOp*
_output_shapes
: *
dtype0*
shape
:
*1
shared_name" Adam_1/discriminator_l1/kernel/m

4Adam_1/discriminator_l1/kernel/m/Read/ReadVariableOpReadVariableOp Adam_1/discriminator_l1/kernel/m*
_output_shapes

:
*
dtype0

Adam_1/discriminator_l1/bias/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:
*/
shared_name Adam_1/discriminator_l1/bias/m

2Adam_1/discriminator_l1/bias/m/Read/ReadVariableOpReadVariableOpAdam_1/discriminator_l1/bias/m*
_output_shapes
:
*
dtype0

 Adam_1/discriminator_l2/kernel/mVarHandleOp*
_output_shapes
: *
dtype0*
shape
:
*1
shared_name" Adam_1/discriminator_l2/kernel/m

4Adam_1/discriminator_l2/kernel/m/Read/ReadVariableOpReadVariableOp Adam_1/discriminator_l2/kernel/m*
_output_shapes

:
*
dtype0

Adam_1/discriminator_l2/bias/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:*/
shared_name Adam_1/discriminator_l2/bias/m

2Adam_1/discriminator_l2/bias/m/Read/ReadVariableOpReadVariableOpAdam_1/discriminator_l2/bias/m*
_output_shapes
:*
dtype0

 Adam_1/discriminator_l3/kernel/mVarHandleOp*
_output_shapes
: *
dtype0*
shape
:*1
shared_name" Adam_1/discriminator_l3/kernel/m

4Adam_1/discriminator_l3/kernel/m/Read/ReadVariableOpReadVariableOp Adam_1/discriminator_l3/kernel/m*
_output_shapes

:*
dtype0

Adam_1/discriminator_l3/bias/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:*/
shared_name Adam_1/discriminator_l3/bias/m

2Adam_1/discriminator_l3/bias/m/Read/ReadVariableOpReadVariableOpAdam_1/discriminator_l3/bias/m*
_output_shapes
:*
dtype0
Є
$Adam_1/discriminator_output/kernel/mVarHandleOp*
_output_shapes
: *
dtype0*
shape
:*5
shared_name&$Adam_1/discriminator_output/kernel/m

8Adam_1/discriminator_output/kernel/m/Read/ReadVariableOpReadVariableOp$Adam_1/discriminator_output/kernel/m*
_output_shapes

:*
dtype0

"Adam_1/discriminator_output/bias/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:*3
shared_name$"Adam_1/discriminator_output/bias/m

6Adam_1/discriminator_output/bias/m/Read/ReadVariableOpReadVariableOp"Adam_1/discriminator_output/bias/m*
_output_shapes
:*
dtype0

Adam/generator_l1/kernel/vVarHandleOp*
_output_shapes
: *
dtype0*
shape
:
*+
shared_nameAdam/generator_l1/kernel/v

.Adam/generator_l1/kernel/v/Read/ReadVariableOpReadVariableOpAdam/generator_l1/kernel/v*
_output_shapes

:
*
dtype0

Adam/generator_l1/bias/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:
*)
shared_nameAdam/generator_l1/bias/v

,Adam/generator_l1/bias/v/Read/ReadVariableOpReadVariableOpAdam/generator_l1/bias/v*
_output_shapes
:
*
dtype0

Adam/generator_l2/kernel/vVarHandleOp*
_output_shapes
: *
dtype0*
shape
:
*+
shared_nameAdam/generator_l2/kernel/v

.Adam/generator_l2/kernel/v/Read/ReadVariableOpReadVariableOpAdam/generator_l2/kernel/v*
_output_shapes

:
*
dtype0

Adam/generator_l2/bias/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:*)
shared_nameAdam/generator_l2/bias/v

,Adam/generator_l2/bias/v/Read/ReadVariableOpReadVariableOpAdam/generator_l2/bias/v*
_output_shapes
:*
dtype0

Adam/generator_l3/kernel/vVarHandleOp*
_output_shapes
: *
dtype0*
shape
:*+
shared_nameAdam/generator_l3/kernel/v

.Adam/generator_l3/kernel/v/Read/ReadVariableOpReadVariableOpAdam/generator_l3/kernel/v*
_output_shapes

:*
dtype0

Adam/generator_l3/bias/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:*)
shared_nameAdam/generator_l3/bias/v

,Adam/generator_l3/bias/v/Read/ReadVariableOpReadVariableOpAdam/generator_l3/bias/v*
_output_shapes
:*
dtype0

Adam/generator_output/kernel/vVarHandleOp*
_output_shapes
: *
dtype0*
shape
:*/
shared_name Adam/generator_output/kernel/v

2Adam/generator_output/kernel/v/Read/ReadVariableOpReadVariableOpAdam/generator_output/kernel/v*
_output_shapes

:*
dtype0

Adam/generator_output/bias/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:*-
shared_nameAdam/generator_output/bias/v

0Adam/generator_output/bias/v/Read/ReadVariableOpReadVariableOpAdam/generator_output/bias/v*
_output_shapes
:*
dtype0

 Adam_1/discriminator_l1/kernel/vVarHandleOp*
_output_shapes
: *
dtype0*
shape
:
*1
shared_name" Adam_1/discriminator_l1/kernel/v

4Adam_1/discriminator_l1/kernel/v/Read/ReadVariableOpReadVariableOp Adam_1/discriminator_l1/kernel/v*
_output_shapes

:
*
dtype0

Adam_1/discriminator_l1/bias/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:
*/
shared_name Adam_1/discriminator_l1/bias/v

2Adam_1/discriminator_l1/bias/v/Read/ReadVariableOpReadVariableOpAdam_1/discriminator_l1/bias/v*
_output_shapes
:
*
dtype0

 Adam_1/discriminator_l2/kernel/vVarHandleOp*
_output_shapes
: *
dtype0*
shape
:
*1
shared_name" Adam_1/discriminator_l2/kernel/v

4Adam_1/discriminator_l2/kernel/v/Read/ReadVariableOpReadVariableOp Adam_1/discriminator_l2/kernel/v*
_output_shapes

:
*
dtype0

Adam_1/discriminator_l2/bias/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:*/
shared_name Adam_1/discriminator_l2/bias/v

2Adam_1/discriminator_l2/bias/v/Read/ReadVariableOpReadVariableOpAdam_1/discriminator_l2/bias/v*
_output_shapes
:*
dtype0

 Adam_1/discriminator_l3/kernel/vVarHandleOp*
_output_shapes
: *
dtype0*
shape
:*1
shared_name" Adam_1/discriminator_l3/kernel/v

4Adam_1/discriminator_l3/kernel/v/Read/ReadVariableOpReadVariableOp Adam_1/discriminator_l3/kernel/v*
_output_shapes

:*
dtype0

Adam_1/discriminator_l3/bias/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:*/
shared_name Adam_1/discriminator_l3/bias/v

2Adam_1/discriminator_l3/bias/v/Read/ReadVariableOpReadVariableOpAdam_1/discriminator_l3/bias/v*
_output_shapes
:*
dtype0
Є
$Adam_1/discriminator_output/kernel/vVarHandleOp*
_output_shapes
: *
dtype0*
shape
:*5
shared_name&$Adam_1/discriminator_output/kernel/v

8Adam_1/discriminator_output/kernel/v/Read/ReadVariableOpReadVariableOp$Adam_1/discriminator_output/kernel/v*
_output_shapes

:*
dtype0

"Adam_1/discriminator_output/bias/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:*3
shared_name$"Adam_1/discriminator_output/bias/v

6Adam_1/discriminator_output/bias/v/Read/ReadVariableOpReadVariableOp"Adam_1/discriminator_output/bias/v*
_output_shapes
:*
dtype0

NoOpNoOp
X
ConstConst"/device:CPU:0*
_output_shapes
: *
dtype0*НW
valueГWBАW BЉW
г
	optimizer
	generator
discriminator
generator_optimizer
discriminator_optimizer
loss
regularization_losses
trainable_variables
	variables
	keras_api
	
signatures


iter

beta_1

beta_2
	decay
learning_rate!m"m#m$m%m&m'm(m)m*m+m,m-m.m/m0m!v"v#v$v %vЁ&vЂ'vЃ(vЄ)vЅ*vІ+vЇ,vЈ-vЉ.vЊ/vЋ0vЌ
ћ
layer-0
layer_with_weights-0
layer-1
layer_with_weights-1
layer-2
layer_with_weights-2
layer-3
layer_with_weights-3
layer-4
regularization_losses
trainable_variables
	variables
	keras_api
ћ
layer-0
layer_with_weights-0
layer-1
layer_with_weights-1
layer-2
layer_with_weights-2
layer-3
layer_with_weights-3
layer-4
regularization_losses
trainable_variables
	variables
 	keras_api
 
 
v
!0
"1
#2
$3
%4
&5
'6
(7
)8
*9
+10
,11
-12
.13
/14
015
v
!0
"1
#2
$3
%4
&5
'6
(7
)8
*9
+10
,11
-12
.13
/14
015
­
1layer_regularization_losses

2layers
3layer_metrics
4non_trainable_variables
5metrics
regularization_losses
trainable_variables
	variables
 
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
 
h

!kernel
"bias
6regularization_losses
7trainable_variables
8	variables
9	keras_api
h

#kernel
$bias
:regularization_losses
;trainable_variables
<	variables
=	keras_api
h

%kernel
&bias
>regularization_losses
?trainable_variables
@	variables
A	keras_api
h

'kernel
(bias
Bregularization_losses
Ctrainable_variables
D	variables
E	keras_api
 
8
!0
"1
#2
$3
%4
&5
'6
(7
8
!0
"1
#2
$3
%4
&5
'6
(7
­
Flayer_regularization_losses

Glayers
Hlayer_metrics
Inon_trainable_variables
Jmetrics
regularization_losses
trainable_variables
	variables
 
h

)kernel
*bias
Kregularization_losses
Ltrainable_variables
M	variables
N	keras_api
h

+kernel
,bias
Oregularization_losses
Ptrainable_variables
Q	variables
R	keras_api
h

-kernel
.bias
Sregularization_losses
Ttrainable_variables
U	variables
V	keras_api
h

/kernel
0bias
Wregularization_losses
Xtrainable_variables
Y	variables
Z	keras_api
 
8
)0
*1
+2
,3
-4
.5
/6
07
8
)0
*1
+2
,3
-4
.5
/6
07
­
[layer_regularization_losses

\layers
]layer_metrics
^non_trainable_variables
_metrics
regularization_losses
trainable_variables
	variables
YW
VARIABLE_VALUEgenerator_l1/kernel0trainable_variables/0/.ATTRIBUTES/VARIABLE_VALUE
WU
VARIABLE_VALUEgenerator_l1/bias0trainable_variables/1/.ATTRIBUTES/VARIABLE_VALUE
YW
VARIABLE_VALUEgenerator_l2/kernel0trainable_variables/2/.ATTRIBUTES/VARIABLE_VALUE
WU
VARIABLE_VALUEgenerator_l2/bias0trainable_variables/3/.ATTRIBUTES/VARIABLE_VALUE
YW
VARIABLE_VALUEgenerator_l3/kernel0trainable_variables/4/.ATTRIBUTES/VARIABLE_VALUE
WU
VARIABLE_VALUEgenerator_l3/bias0trainable_variables/5/.ATTRIBUTES/VARIABLE_VALUE
][
VARIABLE_VALUEgenerator_output/kernel0trainable_variables/6/.ATTRIBUTES/VARIABLE_VALUE
[Y
VARIABLE_VALUEgenerator_output/bias0trainable_variables/7/.ATTRIBUTES/VARIABLE_VALUE
][
VARIABLE_VALUEdiscriminator_l1/kernel0trainable_variables/8/.ATTRIBUTES/VARIABLE_VALUE
[Y
VARIABLE_VALUEdiscriminator_l1/bias0trainable_variables/9/.ATTRIBUTES/VARIABLE_VALUE
^\
VARIABLE_VALUEdiscriminator_l2/kernel1trainable_variables/10/.ATTRIBUTES/VARIABLE_VALUE
\Z
VARIABLE_VALUEdiscriminator_l2/bias1trainable_variables/11/.ATTRIBUTES/VARIABLE_VALUE
^\
VARIABLE_VALUEdiscriminator_l3/kernel1trainable_variables/12/.ATTRIBUTES/VARIABLE_VALUE
\Z
VARIABLE_VALUEdiscriminator_l3/bias1trainable_variables/13/.ATTRIBUTES/VARIABLE_VALUE
b`
VARIABLE_VALUEdiscriminator_output/kernel1trainable_variables/14/.ATTRIBUTES/VARIABLE_VALUE
`^
VARIABLE_VALUEdiscriminator_output/bias1trainable_variables/15/.ATTRIBUTES/VARIABLE_VALUE
 

0
1
 
 

`0
 

!0
"1

!0
"1
­

alayers
blayer_regularization_losses
clayer_metrics
dnon_trainable_variables
emetrics
6regularization_losses
7trainable_variables
8	variables
 

#0
$1

#0
$1
­

flayers
glayer_regularization_losses
hlayer_metrics
inon_trainable_variables
jmetrics
:regularization_losses
;trainable_variables
<	variables
 

%0
&1

%0
&1
­

klayers
llayer_regularization_losses
mlayer_metrics
nnon_trainable_variables
ometrics
>regularization_losses
?trainable_variables
@	variables
 

'0
(1

'0
(1
­

players
qlayer_regularization_losses
rlayer_metrics
snon_trainable_variables
tmetrics
Bregularization_losses
Ctrainable_variables
D	variables
 
#
0
1
2
3
4
 
 
 
 

)0
*1

)0
*1
­

ulayers
vlayer_regularization_losses
wlayer_metrics
xnon_trainable_variables
ymetrics
Kregularization_losses
Ltrainable_variables
M	variables
 

+0
,1

+0
,1
­

zlayers
{layer_regularization_losses
|layer_metrics
}non_trainable_variables
~metrics
Oregularization_losses
Ptrainable_variables
Q	variables
 

-0
.1

-0
.1
Б

layers
 layer_regularization_losses
layer_metrics
non_trainable_variables
metrics
Sregularization_losses
Ttrainable_variables
U	variables
 

/0
01

/0
01
В
layers
 layer_regularization_losses
layer_metrics
non_trainable_variables
metrics
Wregularization_losses
Xtrainable_variables
Y	variables
 
#
0
1
2
3
4
 
 
 
8

total

count
	variables
	keras_api
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
0
1

	variables
|z
VARIABLE_VALUEAdam/generator_l1/kernel/mLtrainable_variables/0/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
zx
VARIABLE_VALUEAdam/generator_l1/bias/mLtrainable_variables/1/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
|z
VARIABLE_VALUEAdam/generator_l2/kernel/mLtrainable_variables/2/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
zx
VARIABLE_VALUEAdam/generator_l2/bias/mLtrainable_variables/3/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
|z
VARIABLE_VALUEAdam/generator_l3/kernel/mLtrainable_variables/4/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
zx
VARIABLE_VALUEAdam/generator_l3/bias/mLtrainable_variables/5/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
~
VARIABLE_VALUEAdam/generator_output/kernel/mLtrainable_variables/6/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
~|
VARIABLE_VALUEAdam/generator_output/bias/mLtrainable_variables/7/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE

VARIABLE_VALUE Adam_1/discriminator_l1/kernel/mLtrainable_variables/8/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
~
VARIABLE_VALUEAdam_1/discriminator_l1/bias/mLtrainable_variables/9/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE

VARIABLE_VALUE Adam_1/discriminator_l2/kernel/mMtrainable_variables/10/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE

VARIABLE_VALUEAdam_1/discriminator_l2/bias/mMtrainable_variables/11/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE

VARIABLE_VALUE Adam_1/discriminator_l3/kernel/mMtrainable_variables/12/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE

VARIABLE_VALUEAdam_1/discriminator_l3/bias/mMtrainable_variables/13/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE

VARIABLE_VALUE$Adam_1/discriminator_output/kernel/mMtrainable_variables/14/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE

VARIABLE_VALUE"Adam_1/discriminator_output/bias/mMtrainable_variables/15/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
|z
VARIABLE_VALUEAdam/generator_l1/kernel/vLtrainable_variables/0/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
zx
VARIABLE_VALUEAdam/generator_l1/bias/vLtrainable_variables/1/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
|z
VARIABLE_VALUEAdam/generator_l2/kernel/vLtrainable_variables/2/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
zx
VARIABLE_VALUEAdam/generator_l2/bias/vLtrainable_variables/3/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
|z
VARIABLE_VALUEAdam/generator_l3/kernel/vLtrainable_variables/4/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
zx
VARIABLE_VALUEAdam/generator_l3/bias/vLtrainable_variables/5/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
~
VARIABLE_VALUEAdam/generator_output/kernel/vLtrainable_variables/6/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
~|
VARIABLE_VALUEAdam/generator_output/bias/vLtrainable_variables/7/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE

VARIABLE_VALUE Adam_1/discriminator_l1/kernel/vLtrainable_variables/8/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
~
VARIABLE_VALUEAdam_1/discriminator_l1/bias/vLtrainable_variables/9/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE

VARIABLE_VALUE Adam_1/discriminator_l2/kernel/vMtrainable_variables/10/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE

VARIABLE_VALUEAdam_1/discriminator_l2/bias/vMtrainable_variables/11/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE

VARIABLE_VALUE Adam_1/discriminator_l3/kernel/vMtrainable_variables/12/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE

VARIABLE_VALUEAdam_1/discriminator_l3/bias/vMtrainable_variables/13/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE

VARIABLE_VALUE$Adam_1/discriminator_output/kernel/vMtrainable_variables/14/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE

VARIABLE_VALUE"Adam_1/discriminator_output/bias/vMtrainable_variables/15/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
z
serving_default_input_1Placeholder*'
_output_shapes
:џџџџџџџџџ*
dtype0*
shape:џџџџџџџџџ

StatefulPartitionedCallStatefulPartitionedCallserving_default_input_1discriminator_l1/kerneldiscriminator_l1/biasdiscriminator_l2/kerneldiscriminator_l2/biasdiscriminator_l3/kerneldiscriminator_l3/biasdiscriminator_output/kerneldiscriminator_output/bias*
Tin
2	*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:џџџџџџџџџ**
_read_only_resource_inputs

*-
config_proto

CPU

GPU 2J 8 *,
f'R%
#__inference_signature_wrapper_71381
O
saver_filenamePlaceholder*
_output_shapes
: *
dtype0*
shape: 
Р
StatefulPartitionedCall_1StatefulPartitionedCallsaver_filenameAdam/iter/Read/ReadVariableOpAdam/beta_1/Read/ReadVariableOpAdam/beta_2/Read/ReadVariableOpAdam/decay/Read/ReadVariableOp&Adam/learning_rate/Read/ReadVariableOp'generator_l1/kernel/Read/ReadVariableOp%generator_l1/bias/Read/ReadVariableOp'generator_l2/kernel/Read/ReadVariableOp%generator_l2/bias/Read/ReadVariableOp'generator_l3/kernel/Read/ReadVariableOp%generator_l3/bias/Read/ReadVariableOp+generator_output/kernel/Read/ReadVariableOp)generator_output/bias/Read/ReadVariableOp+discriminator_l1/kernel/Read/ReadVariableOp)discriminator_l1/bias/Read/ReadVariableOp+discriminator_l2/kernel/Read/ReadVariableOp)discriminator_l2/bias/Read/ReadVariableOp+discriminator_l3/kernel/Read/ReadVariableOp)discriminator_l3/bias/Read/ReadVariableOp/discriminator_output/kernel/Read/ReadVariableOp-discriminator_output/bias/Read/ReadVariableOptotal/Read/ReadVariableOpcount/Read/ReadVariableOp.Adam/generator_l1/kernel/m/Read/ReadVariableOp,Adam/generator_l1/bias/m/Read/ReadVariableOp.Adam/generator_l2/kernel/m/Read/ReadVariableOp,Adam/generator_l2/bias/m/Read/ReadVariableOp.Adam/generator_l3/kernel/m/Read/ReadVariableOp,Adam/generator_l3/bias/m/Read/ReadVariableOp2Adam/generator_output/kernel/m/Read/ReadVariableOp0Adam/generator_output/bias/m/Read/ReadVariableOp4Adam_1/discriminator_l1/kernel/m/Read/ReadVariableOp2Adam_1/discriminator_l1/bias/m/Read/ReadVariableOp4Adam_1/discriminator_l2/kernel/m/Read/ReadVariableOp2Adam_1/discriminator_l2/bias/m/Read/ReadVariableOp4Adam_1/discriminator_l3/kernel/m/Read/ReadVariableOp2Adam_1/discriminator_l3/bias/m/Read/ReadVariableOp8Adam_1/discriminator_output/kernel/m/Read/ReadVariableOp6Adam_1/discriminator_output/bias/m/Read/ReadVariableOp.Adam/generator_l1/kernel/v/Read/ReadVariableOp,Adam/generator_l1/bias/v/Read/ReadVariableOp.Adam/generator_l2/kernel/v/Read/ReadVariableOp,Adam/generator_l2/bias/v/Read/ReadVariableOp.Adam/generator_l3/kernel/v/Read/ReadVariableOp,Adam/generator_l3/bias/v/Read/ReadVariableOp2Adam/generator_output/kernel/v/Read/ReadVariableOp0Adam/generator_output/bias/v/Read/ReadVariableOp4Adam_1/discriminator_l1/kernel/v/Read/ReadVariableOp2Adam_1/discriminator_l1/bias/v/Read/ReadVariableOp4Adam_1/discriminator_l2/kernel/v/Read/ReadVariableOp2Adam_1/discriminator_l2/bias/v/Read/ReadVariableOp4Adam_1/discriminator_l3/kernel/v/Read/ReadVariableOp2Adam_1/discriminator_l3/bias/v/Read/ReadVariableOp8Adam_1/discriminator_output/kernel/v/Read/ReadVariableOp6Adam_1/discriminator_output/bias/v/Read/ReadVariableOpConst*D
Tin=
;29	*
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
GPU 2J 8 *'
f"R 
__inference__traced_save_72153
я
StatefulPartitionedCall_2StatefulPartitionedCallsaver_filename	Adam/iterAdam/beta_1Adam/beta_2
Adam/decayAdam/learning_rategenerator_l1/kernelgenerator_l1/biasgenerator_l2/kernelgenerator_l2/biasgenerator_l3/kernelgenerator_l3/biasgenerator_output/kernelgenerator_output/biasdiscriminator_l1/kerneldiscriminator_l1/biasdiscriminator_l2/kerneldiscriminator_l2/biasdiscriminator_l3/kerneldiscriminator_l3/biasdiscriminator_output/kerneldiscriminator_output/biastotalcountAdam/generator_l1/kernel/mAdam/generator_l1/bias/mAdam/generator_l2/kernel/mAdam/generator_l2/bias/mAdam/generator_l3/kernel/mAdam/generator_l3/bias/mAdam/generator_output/kernel/mAdam/generator_output/bias/m Adam_1/discriminator_l1/kernel/mAdam_1/discriminator_l1/bias/m Adam_1/discriminator_l2/kernel/mAdam_1/discriminator_l2/bias/m Adam_1/discriminator_l3/kernel/mAdam_1/discriminator_l3/bias/m$Adam_1/discriminator_output/kernel/m"Adam_1/discriminator_output/bias/mAdam/generator_l1/kernel/vAdam/generator_l1/bias/vAdam/generator_l2/kernel/vAdam/generator_l2/bias/vAdam/generator_l3/kernel/vAdam/generator_l3/bias/vAdam/generator_output/kernel/vAdam/generator_output/bias/v Adam_1/discriminator_l1/kernel/vAdam_1/discriminator_l1/bias/v Adam_1/discriminator_l2/kernel/vAdam_1/discriminator_l2/bias/v Adam_1/discriminator_l3/kernel/vAdam_1/discriminator_l3/bias/v$Adam_1/discriminator_output/kernel/v"Adam_1/discriminator_output/bias/v*C
Tin<
:28*
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
GPU 2J 8 **
f%R#
!__inference__traced_restore_72328Ел

л
џ
H__inference_discriminator_layer_call_and_return_conditional_losses_71116

inputs
discriminator_l1_71095
discriminator_l1_71097
discriminator_l2_71100
discriminator_l2_71102
discriminator_l3_71105
discriminator_l3_71107
discriminator_output_71110
discriminator_output_71112
identityЂ(discriminator_l1/StatefulPartitionedCallЂ(discriminator_l2/StatefulPartitionedCallЂ(discriminator_l3/StatefulPartitionedCallЂ,discriminator_output/StatefulPartitionedCallЙ
(discriminator_l1/StatefulPartitionedCallStatefulPartitionedCallinputsdiscriminator_l1_71095discriminator_l1_71097*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:џџџџџџџџџ
*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8 *T
fORM
K__inference_discriminator_l1_layer_call_and_return_conditional_losses_709672*
(discriminator_l1/StatefulPartitionedCallф
(discriminator_l2/StatefulPartitionedCallStatefulPartitionedCall1discriminator_l1/StatefulPartitionedCall:output:0discriminator_l2_71100discriminator_l2_71102*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:џџџџџџџџџ*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8 *T
fORM
K__inference_discriminator_l2_layer_call_and_return_conditional_losses_709942*
(discriminator_l2/StatefulPartitionedCallф
(discriminator_l3/StatefulPartitionedCallStatefulPartitionedCall1discriminator_l2/StatefulPartitionedCall:output:0discriminator_l3_71105discriminator_l3_71107*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:џџџџџџџџџ*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8 *T
fORM
K__inference_discriminator_l3_layer_call_and_return_conditional_losses_710212*
(discriminator_l3/StatefulPartitionedCallј
,discriminator_output/StatefulPartitionedCallStatefulPartitionedCall1discriminator_l3/StatefulPartitionedCall:output:0discriminator_output_71110discriminator_output_71112*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:џџџџџџџџџ*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8 *X
fSRQ
O__inference_discriminator_output_layer_call_and_return_conditional_losses_710482.
,discriminator_output/StatefulPartitionedCallЙ
IdentityIdentity5discriminator_output/StatefulPartitionedCall:output:0)^discriminator_l1/StatefulPartitionedCall)^discriminator_l2/StatefulPartitionedCall)^discriminator_l3/StatefulPartitionedCall-^discriminator_output/StatefulPartitionedCall*
T0*'
_output_shapes
:џџџџџџџџџ2

Identity"
identityIdentity:output:0*F
_input_shapes5
3:џџџџџџџџџ::::::::2T
(discriminator_l1/StatefulPartitionedCall(discriminator_l1/StatefulPartitionedCall2T
(discriminator_l2/StatefulPartitionedCall(discriminator_l2/StatefulPartitionedCall2T
(discriminator_l3/StatefulPartitionedCall(discriminator_l3/StatefulPartitionedCall2\
,discriminator_output/StatefulPartitionedCall,discriminator_output/StatefulPartitionedCall:O K
'
_output_shapes
:џџџџџџџџџ
 
_user_specified_nameinputs
ћ	
ш
O__inference_discriminator_output_layer_call_and_return_conditional_losses_71956

inputs"
matmul_readvariableop_resource#
biasadd_readvariableop_resource
identityЂBiasAdd/ReadVariableOpЂMatMul/ReadVariableOp
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:*
dtype02
MatMul/ReadVariableOps
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2
MatMul
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:*
dtype02
BiasAdd/ReadVariableOp
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2	
BiasAdda
SigmoidSigmoidBiasAdd:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2	
Sigmoid
IdentityIdentitySigmoid:y:0^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*
T0*'
_output_shapes
:џџџџџџџџџ2

Identity"
identityIdentity:output:0*.
_input_shapes
:џџџџџџџџџ::20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:O K
'
_output_shapes
:џџџџџџџџџ
 
_user_specified_nameinputs
Ь
щ
-__inference_discriminator_layer_call_fn_71135
discriminator_input
unknown
	unknown_0
	unknown_1
	unknown_2
	unknown_3
	unknown_4
	unknown_5
	unknown_6
identityЂStatefulPartitionedCallг
StatefulPartitionedCallStatefulPartitionedCalldiscriminator_inputunknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6*
Tin
2	*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:џџџџџџџџџ**
_read_only_resource_inputs

*-
config_proto

CPU

GPU 2J 8 *Q
fLRJ
H__inference_discriminator_layer_call_and_return_conditional_losses_711162
StatefulPartitionedCall
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*'
_output_shapes
:џџџџџџџџџ2

Identity"
identityIdentity:output:0*F
_input_shapes5
3:џџџџџџџџџ::::::::22
StatefulPartitionedCallStatefulPartitionedCall:\ X
'
_output_shapes
:џџџџџџџџџ
-
_user_specified_namediscriminator_input
у

,__inference_generator_l1_layer_call_fn_71825

inputs
unknown
	unknown_0
identityЂStatefulPartitionedCallї
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:џџџџџџџџџ
*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8 *P
fKRI
G__inference_generator_l1_layer_call_and_return_conditional_losses_707392
StatefulPartitionedCall
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*'
_output_shapes
:џџџџџџџџџ
2

Identity"
identityIdentity:output:0*.
_input_shapes
:џџџџџџџџџ::22
StatefulPartitionedCallStatefulPartitionedCall:O K
'
_output_shapes
:џџџџџџџџџ
 
_user_specified_nameinputs
Ѕ
м
-__inference_discriminator_layer_call_fn_71805

inputs
unknown
	unknown_0
	unknown_1
	unknown_2
	unknown_3
	unknown_4
	unknown_5
	unknown_6
identityЂStatefulPartitionedCallЦ
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6*
Tin
2	*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:џџџџџџџџџ**
_read_only_resource_inputs

*-
config_proto

CPU

GPU 2J 8 *Q
fLRJ
H__inference_discriminator_layer_call_and_return_conditional_losses_711612
StatefulPartitionedCall
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*'
_output_shapes
:џџџџџџџџџ2

Identity"
identityIdentity:output:0*F
_input_shapes5
3:џџџџџџџџџ::::::::22
StatefulPartitionedCallStatefulPartitionedCall:O K
'
_output_shapes
:џџџџџџџџџ
 
_user_specified_nameinputs
ѕ	
ф
K__inference_generator_output_layer_call_and_return_conditional_losses_71876

inputs"
matmul_readvariableop_resource#
biasadd_readvariableop_resource
identityЂBiasAdd/ReadVariableOpЂMatMul/ReadVariableOp
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:*
dtype02
MatMul/ReadVariableOps
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2
MatMul
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:*
dtype02
BiasAdd/ReadVariableOp
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2	
BiasAddX
ReluReluBiasAdd:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
Relu
IdentityIdentityRelu:activations:0^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*
T0*'
_output_shapes
:џџџџџџџџџ2

Identity"
identityIdentity:output:0*.
_input_shapes
:џџџџџџџџџ::20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:O K
'
_output_shapes
:џџџџџџџџџ
 
_user_specified_nameinputs
И
с
)__inference_generator_layer_call_fn_70952
generator_input
unknown
	unknown_0
	unknown_1
	unknown_2
	unknown_3
	unknown_4
	unknown_5
	unknown_6
identityЂStatefulPartitionedCallЫ
StatefulPartitionedCallStatefulPartitionedCallgenerator_inputunknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6*
Tin
2	*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:џџџџџџџџџ**
_read_only_resource_inputs

*-
config_proto

CPU

GPU 2J 8 *M
fHRF
D__inference_generator_layer_call_and_return_conditional_losses_709332
StatefulPartitionedCall
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*'
_output_shapes
:џџџџџџџџџ2

Identity"
identityIdentity:output:0*F
_input_shapes5
3:џџџџџџџџџ::::::::22
StatefulPartitionedCallStatefulPartitionedCall:X T
'
_output_shapes
:џџџџџџџџџ
)
_user_specified_namegenerator_input
у

,__inference_generator_l3_layer_call_fn_71865

inputs
unknown
	unknown_0
identityЂStatefulPartitionedCallї
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:џџџџџџџџџ*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8 *P
fKRI
G__inference_generator_l3_layer_call_and_return_conditional_losses_707932
StatefulPartitionedCall
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*'
_output_shapes
:џџџџџџџџџ2

Identity"
identityIdentity:output:0*.
_input_shapes
:џџџџџџџџџ::22
StatefulPartitionedCallStatefulPartitionedCall:O K
'
_output_shapes
:џџџџџџџџџ
 
_user_specified_nameinputs
В
д
D__inference_generator_layer_call_and_return_conditional_losses_70837
generator_input
generator_l1_70750
generator_l1_70752
generator_l2_70777
generator_l2_70779
generator_l3_70804
generator_l3_70806
generator_output_70831
generator_output_70833
identityЂ$generator_l1/StatefulPartitionedCallЂ$generator_l2/StatefulPartitionedCallЂ$generator_l3/StatefulPartitionedCallЂ(generator_output/StatefulPartitionedCallЎ
$generator_l1/StatefulPartitionedCallStatefulPartitionedCallgenerator_inputgenerator_l1_70750generator_l1_70752*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:џџџџџџџџџ
*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8 *P
fKRI
G__inference_generator_l1_layer_call_and_return_conditional_losses_707392&
$generator_l1/StatefulPartitionedCallЬ
$generator_l2/StatefulPartitionedCallStatefulPartitionedCall-generator_l1/StatefulPartitionedCall:output:0generator_l2_70777generator_l2_70779*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:џџџџџџџџџ*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8 *P
fKRI
G__inference_generator_l2_layer_call_and_return_conditional_losses_707662&
$generator_l2/StatefulPartitionedCallЬ
$generator_l3/StatefulPartitionedCallStatefulPartitionedCall-generator_l2/StatefulPartitionedCall:output:0generator_l3_70804generator_l3_70806*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:џџџџџџџџџ*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8 *P
fKRI
G__inference_generator_l3_layer_call_and_return_conditional_losses_707932&
$generator_l3/StatefulPartitionedCallр
(generator_output/StatefulPartitionedCallStatefulPartitionedCall-generator_l3/StatefulPartitionedCall:output:0generator_output_70831generator_output_70833*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:џџџџџџџџџ*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8 *T
fORM
K__inference_generator_output_layer_call_and_return_conditional_losses_708202*
(generator_output/StatefulPartitionedCallЅ
IdentityIdentity1generator_output/StatefulPartitionedCall:output:0%^generator_l1/StatefulPartitionedCall%^generator_l2/StatefulPartitionedCall%^generator_l3/StatefulPartitionedCall)^generator_output/StatefulPartitionedCall*
T0*'
_output_shapes
:џџџџџџџџџ2

Identity"
identityIdentity:output:0*F
_input_shapes5
3:џџџџџџџџџ::::::::2L
$generator_l1/StatefulPartitionedCall$generator_l1/StatefulPartitionedCall2L
$generator_l2/StatefulPartitionedCall$generator_l2/StatefulPartitionedCall2L
$generator_l3/StatefulPartitionedCall$generator_l3/StatefulPartitionedCall2T
(generator_output/StatefulPartitionedCall(generator_output/StatefulPartitionedCall:X T
'
_output_shapes
:џџџџџџџџџ
)
_user_specified_namegenerator_input
ё	
р
G__inference_generator_l1_layer_call_and_return_conditional_losses_71816

inputs"
matmul_readvariableop_resource#
biasadd_readvariableop_resource
identityЂBiasAdd/ReadVariableOpЂMatMul/ReadVariableOp
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:
*
dtype02
MatMul/ReadVariableOps
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ
2
MatMul
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:
*
dtype02
BiasAdd/ReadVariableOp
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ
2	
BiasAddX
ReluReluBiasAdd:output:0*
T0*'
_output_shapes
:џџџџџџџџџ
2
Relu
IdentityIdentityRelu:activations:0^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*
T0*'
_output_shapes
:џџџџџџџџџ
2

Identity"
identityIdentity:output:0*.
_input_shapes
:џџџџџџџџџ::20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:O K
'
_output_shapes
:џџџџџџџџџ
 
_user_specified_nameinputs
ѕ	
ф
K__inference_discriminator_l1_layer_call_and_return_conditional_losses_71896

inputs"
matmul_readvariableop_resource#
biasadd_readvariableop_resource
identityЂBiasAdd/ReadVariableOpЂMatMul/ReadVariableOp
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:
*
dtype02
MatMul/ReadVariableOps
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ
2
MatMul
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:
*
dtype02
BiasAdd/ReadVariableOp
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ
2	
BiasAddX
ReluReluBiasAdd:output:0*
T0*'
_output_shapes
:џџџџџџџџџ
2
Relu
IdentityIdentityRelu:activations:0^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*
T0*'
_output_shapes
:џџџџџџџџџ
2

Identity"
identityIdentity:output:0*.
_input_shapes
:џџџџџџџџџ::20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:O K
'
_output_shapes
:џџџџџџџџџ
 
_user_specified_nameinputs
ё	
р
G__inference_generator_l2_layer_call_and_return_conditional_losses_70766

inputs"
matmul_readvariableop_resource#
biasadd_readvariableop_resource
identityЂBiasAdd/ReadVariableOpЂMatMul/ReadVariableOp
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:
*
dtype02
MatMul/ReadVariableOps
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2
MatMul
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:*
dtype02
BiasAdd/ReadVariableOp
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2	
BiasAddX
ReluReluBiasAdd:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
Relu
IdentityIdentityRelu:activations:0^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*
T0*'
_output_shapes
:џџџџџџџџџ2

Identity"
identityIdentity:output:0*.
_input_shapes
:џџџџџџџџџ
::20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:O K
'
_output_shapes
:џџџџџџџџџ

 
_user_specified_nameinputs
ђ@
Ј	
 __inference__wrapped_model_70724
input_1I
Ecgan_31_discriminator_discriminator_l1_matmul_readvariableop_resourceJ
Fcgan_31_discriminator_discriminator_l1_biasadd_readvariableop_resourceI
Ecgan_31_discriminator_discriminator_l2_matmul_readvariableop_resourceJ
Fcgan_31_discriminator_discriminator_l2_biasadd_readvariableop_resourceI
Ecgan_31_discriminator_discriminator_l3_matmul_readvariableop_resourceJ
Fcgan_31_discriminator_discriminator_l3_biasadd_readvariableop_resourceM
Icgan_31_discriminator_discriminator_output_matmul_readvariableop_resourceN
Jcgan_31_discriminator_discriminator_output_biasadd_readvariableop_resource
identityЂ=cgan_31/discriminator/discriminator_l1/BiasAdd/ReadVariableOpЂ<cgan_31/discriminator/discriminator_l1/MatMul/ReadVariableOpЂ=cgan_31/discriminator/discriminator_l2/BiasAdd/ReadVariableOpЂ<cgan_31/discriminator/discriminator_l2/MatMul/ReadVariableOpЂ=cgan_31/discriminator/discriminator_l3/BiasAdd/ReadVariableOpЂ<cgan_31/discriminator/discriminator_l3/MatMul/ReadVariableOpЂAcgan_31/discriminator/discriminator_output/BiasAdd/ReadVariableOpЂ@cgan_31/discriminator/discriminator_output/MatMul/ReadVariableOp
<cgan_31/discriminator/discriminator_l1/MatMul/ReadVariableOpReadVariableOpEcgan_31_discriminator_discriminator_l1_matmul_readvariableop_resource*
_output_shapes

:
*
dtype02>
<cgan_31/discriminator/discriminator_l1/MatMul/ReadVariableOpщ
-cgan_31/discriminator/discriminator_l1/MatMulMatMulinput_1Dcgan_31/discriminator/discriminator_l1/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ
2/
-cgan_31/discriminator/discriminator_l1/MatMul
=cgan_31/discriminator/discriminator_l1/BiasAdd/ReadVariableOpReadVariableOpFcgan_31_discriminator_discriminator_l1_biasadd_readvariableop_resource*
_output_shapes
:
*
dtype02?
=cgan_31/discriminator/discriminator_l1/BiasAdd/ReadVariableOp
.cgan_31/discriminator/discriminator_l1/BiasAddBiasAdd7cgan_31/discriminator/discriminator_l1/MatMul:product:0Ecgan_31/discriminator/discriminator_l1/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ
20
.cgan_31/discriminator/discriminator_l1/BiasAddЭ
+cgan_31/discriminator/discriminator_l1/ReluRelu7cgan_31/discriminator/discriminator_l1/BiasAdd:output:0*
T0*'
_output_shapes
:џџџџџџџџџ
2-
+cgan_31/discriminator/discriminator_l1/Relu
<cgan_31/discriminator/discriminator_l2/MatMul/ReadVariableOpReadVariableOpEcgan_31_discriminator_discriminator_l2_matmul_readvariableop_resource*
_output_shapes

:
*
dtype02>
<cgan_31/discriminator/discriminator_l2/MatMul/ReadVariableOp
-cgan_31/discriminator/discriminator_l2/MatMulMatMul9cgan_31/discriminator/discriminator_l1/Relu:activations:0Dcgan_31/discriminator/discriminator_l2/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2/
-cgan_31/discriminator/discriminator_l2/MatMul
=cgan_31/discriminator/discriminator_l2/BiasAdd/ReadVariableOpReadVariableOpFcgan_31_discriminator_discriminator_l2_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02?
=cgan_31/discriminator/discriminator_l2/BiasAdd/ReadVariableOp
.cgan_31/discriminator/discriminator_l2/BiasAddBiasAdd7cgan_31/discriminator/discriminator_l2/MatMul:product:0Ecgan_31/discriminator/discriminator_l2/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ20
.cgan_31/discriminator/discriminator_l2/BiasAddЭ
+cgan_31/discriminator/discriminator_l2/ReluRelu7cgan_31/discriminator/discriminator_l2/BiasAdd:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2-
+cgan_31/discriminator/discriminator_l2/Relu
<cgan_31/discriminator/discriminator_l3/MatMul/ReadVariableOpReadVariableOpEcgan_31_discriminator_discriminator_l3_matmul_readvariableop_resource*
_output_shapes

:*
dtype02>
<cgan_31/discriminator/discriminator_l3/MatMul/ReadVariableOp
-cgan_31/discriminator/discriminator_l3/MatMulMatMul9cgan_31/discriminator/discriminator_l2/Relu:activations:0Dcgan_31/discriminator/discriminator_l3/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2/
-cgan_31/discriminator/discriminator_l3/MatMul
=cgan_31/discriminator/discriminator_l3/BiasAdd/ReadVariableOpReadVariableOpFcgan_31_discriminator_discriminator_l3_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02?
=cgan_31/discriminator/discriminator_l3/BiasAdd/ReadVariableOp
.cgan_31/discriminator/discriminator_l3/BiasAddBiasAdd7cgan_31/discriminator/discriminator_l3/MatMul:product:0Ecgan_31/discriminator/discriminator_l3/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ20
.cgan_31/discriminator/discriminator_l3/BiasAddЭ
+cgan_31/discriminator/discriminator_l3/ReluRelu7cgan_31/discriminator/discriminator_l3/BiasAdd:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2-
+cgan_31/discriminator/discriminator_l3/Relu
@cgan_31/discriminator/discriminator_output/MatMul/ReadVariableOpReadVariableOpIcgan_31_discriminator_discriminator_output_matmul_readvariableop_resource*
_output_shapes

:*
dtype02B
@cgan_31/discriminator/discriminator_output/MatMul/ReadVariableOpЇ
1cgan_31/discriminator/discriminator_output/MatMulMatMul9cgan_31/discriminator/discriminator_l3/Relu:activations:0Hcgan_31/discriminator/discriminator_output/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ23
1cgan_31/discriminator/discriminator_output/MatMul
Acgan_31/discriminator/discriminator_output/BiasAdd/ReadVariableOpReadVariableOpJcgan_31_discriminator_discriminator_output_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02C
Acgan_31/discriminator/discriminator_output/BiasAdd/ReadVariableOp­
2cgan_31/discriminator/discriminator_output/BiasAddBiasAdd;cgan_31/discriminator/discriminator_output/MatMul:product:0Icgan_31/discriminator/discriminator_output/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ24
2cgan_31/discriminator/discriminator_output/BiasAddт
2cgan_31/discriminator/discriminator_output/SigmoidSigmoid;cgan_31/discriminator/discriminator_output/BiasAdd:output:0*
T0*'
_output_shapes
:џџџџџџџџџ24
2cgan_31/discriminator/discriminator_output/Sigmoid
IdentityIdentity6cgan_31/discriminator/discriminator_output/Sigmoid:y:0>^cgan_31/discriminator/discriminator_l1/BiasAdd/ReadVariableOp=^cgan_31/discriminator/discriminator_l1/MatMul/ReadVariableOp>^cgan_31/discriminator/discriminator_l2/BiasAdd/ReadVariableOp=^cgan_31/discriminator/discriminator_l2/MatMul/ReadVariableOp>^cgan_31/discriminator/discriminator_l3/BiasAdd/ReadVariableOp=^cgan_31/discriminator/discriminator_l3/MatMul/ReadVariableOpB^cgan_31/discriminator/discriminator_output/BiasAdd/ReadVariableOpA^cgan_31/discriminator/discriminator_output/MatMul/ReadVariableOp*
T0*'
_output_shapes
:џџџџџџџџџ2

Identity"
identityIdentity:output:0*F
_input_shapes5
3:џџџџџџџџџ::::::::2~
=cgan_31/discriminator/discriminator_l1/BiasAdd/ReadVariableOp=cgan_31/discriminator/discriminator_l1/BiasAdd/ReadVariableOp2|
<cgan_31/discriminator/discriminator_l1/MatMul/ReadVariableOp<cgan_31/discriminator/discriminator_l1/MatMul/ReadVariableOp2~
=cgan_31/discriminator/discriminator_l2/BiasAdd/ReadVariableOp=cgan_31/discriminator/discriminator_l2/BiasAdd/ReadVariableOp2|
<cgan_31/discriminator/discriminator_l2/MatMul/ReadVariableOp<cgan_31/discriminator/discriminator_l2/MatMul/ReadVariableOp2~
=cgan_31/discriminator/discriminator_l3/BiasAdd/ReadVariableOp=cgan_31/discriminator/discriminator_l3/BiasAdd/ReadVariableOp2|
<cgan_31/discriminator/discriminator_l3/MatMul/ReadVariableOp<cgan_31/discriminator/discriminator_l3/MatMul/ReadVariableOp2
Acgan_31/discriminator/discriminator_output/BiasAdd/ReadVariableOpAcgan_31/discriminator/discriminator_output/BiasAdd/ReadVariableOp2
@cgan_31/discriminator/discriminator_output/MatMul/ReadVariableOp@cgan_31/discriminator/discriminator_output/MatMul/ReadVariableOp:P L
'
_output_shapes
:џџџџџџџџџ
!
_user_specified_name	input_1
И
с
)__inference_generator_layer_call_fn_70907
generator_input
unknown
	unknown_0
	unknown_1
	unknown_2
	unknown_3
	unknown_4
	unknown_5
	unknown_6
identityЂStatefulPartitionedCallЫ
StatefulPartitionedCallStatefulPartitionedCallgenerator_inputunknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6*
Tin
2	*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:џџџџџџџџџ**
_read_only_resource_inputs

*-
config_proto

CPU

GPU 2J 8 *M
fHRF
D__inference_generator_layer_call_and_return_conditional_losses_708882
StatefulPartitionedCall
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*'
_output_shapes
:џџџџџџџџџ2

Identity"
identityIdentity:output:0*F
_input_shapes5
3:џџџџџџџџџ::::::::22
StatefulPartitionedCallStatefulPartitionedCall:X T
'
_output_shapes
:џџџџџџџџџ
)
_user_specified_namegenerator_input
+
Ћ
D__inference_generator_layer_call_and_return_conditional_losses_71657

inputs/
+generator_l1_matmul_readvariableop_resource0
,generator_l1_biasadd_readvariableop_resource/
+generator_l2_matmul_readvariableop_resource0
,generator_l2_biasadd_readvariableop_resource/
+generator_l3_matmul_readvariableop_resource0
,generator_l3_biasadd_readvariableop_resource3
/generator_output_matmul_readvariableop_resource4
0generator_output_biasadd_readvariableop_resource
identityЂ#generator_l1/BiasAdd/ReadVariableOpЂ"generator_l1/MatMul/ReadVariableOpЂ#generator_l2/BiasAdd/ReadVariableOpЂ"generator_l2/MatMul/ReadVariableOpЂ#generator_l3/BiasAdd/ReadVariableOpЂ"generator_l3/MatMul/ReadVariableOpЂ'generator_output/BiasAdd/ReadVariableOpЂ&generator_output/MatMul/ReadVariableOpД
"generator_l1/MatMul/ReadVariableOpReadVariableOp+generator_l1_matmul_readvariableop_resource*
_output_shapes

:
*
dtype02$
"generator_l1/MatMul/ReadVariableOp
generator_l1/MatMulMatMulinputs*generator_l1/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ
2
generator_l1/MatMulГ
#generator_l1/BiasAdd/ReadVariableOpReadVariableOp,generator_l1_biasadd_readvariableop_resource*
_output_shapes
:
*
dtype02%
#generator_l1/BiasAdd/ReadVariableOpЕ
generator_l1/BiasAddBiasAddgenerator_l1/MatMul:product:0+generator_l1/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ
2
generator_l1/BiasAdd
generator_l1/ReluRelugenerator_l1/BiasAdd:output:0*
T0*'
_output_shapes
:џџџџџџџџџ
2
generator_l1/ReluД
"generator_l2/MatMul/ReadVariableOpReadVariableOp+generator_l2_matmul_readvariableop_resource*
_output_shapes

:
*
dtype02$
"generator_l2/MatMul/ReadVariableOpГ
generator_l2/MatMulMatMulgenerator_l1/Relu:activations:0*generator_l2/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2
generator_l2/MatMulГ
#generator_l2/BiasAdd/ReadVariableOpReadVariableOp,generator_l2_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02%
#generator_l2/BiasAdd/ReadVariableOpЕ
generator_l2/BiasAddBiasAddgenerator_l2/MatMul:product:0+generator_l2/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2
generator_l2/BiasAdd
generator_l2/ReluRelugenerator_l2/BiasAdd:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
generator_l2/ReluД
"generator_l3/MatMul/ReadVariableOpReadVariableOp+generator_l3_matmul_readvariableop_resource*
_output_shapes

:*
dtype02$
"generator_l3/MatMul/ReadVariableOpГ
generator_l3/MatMulMatMulgenerator_l2/Relu:activations:0*generator_l3/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2
generator_l3/MatMulГ
#generator_l3/BiasAdd/ReadVariableOpReadVariableOp,generator_l3_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02%
#generator_l3/BiasAdd/ReadVariableOpЕ
generator_l3/BiasAddBiasAddgenerator_l3/MatMul:product:0+generator_l3/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2
generator_l3/BiasAdd
generator_l3/ReluRelugenerator_l3/BiasAdd:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
generator_l3/ReluР
&generator_output/MatMul/ReadVariableOpReadVariableOp/generator_output_matmul_readvariableop_resource*
_output_shapes

:*
dtype02(
&generator_output/MatMul/ReadVariableOpП
generator_output/MatMulMatMulgenerator_l3/Relu:activations:0.generator_output/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2
generator_output/MatMulП
'generator_output/BiasAdd/ReadVariableOpReadVariableOp0generator_output_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02)
'generator_output/BiasAdd/ReadVariableOpХ
generator_output/BiasAddBiasAdd!generator_output/MatMul:product:0/generator_output/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2
generator_output/BiasAdd
generator_output/ReluRelu!generator_output/BiasAdd:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
generator_output/ReluЋ
IdentityIdentity#generator_output/Relu:activations:0$^generator_l1/BiasAdd/ReadVariableOp#^generator_l1/MatMul/ReadVariableOp$^generator_l2/BiasAdd/ReadVariableOp#^generator_l2/MatMul/ReadVariableOp$^generator_l3/BiasAdd/ReadVariableOp#^generator_l3/MatMul/ReadVariableOp(^generator_output/BiasAdd/ReadVariableOp'^generator_output/MatMul/ReadVariableOp*
T0*'
_output_shapes
:џџџџџџџџџ2

Identity"
identityIdentity:output:0*F
_input_shapes5
3:џџџџџџџџџ::::::::2J
#generator_l1/BiasAdd/ReadVariableOp#generator_l1/BiasAdd/ReadVariableOp2H
"generator_l1/MatMul/ReadVariableOp"generator_l1/MatMul/ReadVariableOp2J
#generator_l2/BiasAdd/ReadVariableOp#generator_l2/BiasAdd/ReadVariableOp2H
"generator_l2/MatMul/ReadVariableOp"generator_l2/MatMul/ReadVariableOp2J
#generator_l3/BiasAdd/ReadVariableOp#generator_l3/BiasAdd/ReadVariableOp2H
"generator_l3/MatMul/ReadVariableOp"generator_l3/MatMul/ReadVariableOp2R
'generator_output/BiasAdd/ReadVariableOp'generator_output/BiasAdd/ReadVariableOp2P
&generator_output/MatMul/ReadVariableOp&generator_output/MatMul/ReadVariableOp:O K
'
_output_shapes
:џџџџџџџџџ
 
_user_specified_nameinputs
ѕ	
ф
K__inference_discriminator_l1_layer_call_and_return_conditional_losses_70967

inputs"
matmul_readvariableop_resource#
biasadd_readvariableop_resource
identityЂBiasAdd/ReadVariableOpЂMatMul/ReadVariableOp
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:
*
dtype02
MatMul/ReadVariableOps
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ
2
MatMul
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:
*
dtype02
BiasAdd/ReadVariableOp
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ
2	
BiasAddX
ReluReluBiasAdd:output:0*
T0*'
_output_shapes
:џџџџџџџџџ
2
Relu
IdentityIdentityRelu:activations:0^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*
T0*'
_output_shapes
:џџџџџџџџџ
2

Identity"
identityIdentity:output:0*.
_input_shapes
:џџџџџџџџџ::20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:O K
'
_output_shapes
:џџџџџџџџџ
 
_user_specified_nameinputs
ё	
р
G__inference_generator_l3_layer_call_and_return_conditional_losses_70793

inputs"
matmul_readvariableop_resource#
biasadd_readvariableop_resource
identityЂBiasAdd/ReadVariableOpЂMatMul/ReadVariableOp
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:*
dtype02
MatMul/ReadVariableOps
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2
MatMul
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:*
dtype02
BiasAdd/ReadVariableOp
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2	
BiasAddX
ReluReluBiasAdd:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
Relu
IdentityIdentityRelu:activations:0^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*
T0*'
_output_shapes
:џџџџџџџџџ2

Identity"
identityIdentity:output:0*.
_input_shapes
:џџџџџџџџџ::20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:O K
'
_output_shapes
:џџџџџџџџџ
 
_user_specified_nameinputs
і
г
#__inference_signature_wrapper_71381
input_1
unknown
	unknown_0
	unknown_1
	unknown_2
	unknown_3
	unknown_4
	unknown_5
	unknown_6
identityЂStatefulPartitionedCall
StatefulPartitionedCallStatefulPartitionedCallinput_1unknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6*
Tin
2	*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:џџџџџџџџџ**
_read_only_resource_inputs

*-
config_proto

CPU

GPU 2J 8 *)
f$R"
 __inference__wrapped_model_707242
StatefulPartitionedCall
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*'
_output_shapes
:џџџџџџџџџ2

Identity"
identityIdentity:output:0*F
_input_shapes5
3:џџџџџџџџџ::::::::22
StatefulPartitionedCallStatefulPartitionedCall:P L
'
_output_shapes
:џџџџџџџџџ
!
_user_specified_name	input_1
ё	
р
G__inference_generator_l3_layer_call_and_return_conditional_losses_71856

inputs"
matmul_readvariableop_resource#
biasadd_readvariableop_resource
identityЂBiasAdd/ReadVariableOpЂMatMul/ReadVariableOp
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:*
dtype02
MatMul/ReadVariableOps
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2
MatMul
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:*
dtype02
BiasAdd/ReadVariableOp
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2	
BiasAddX
ReluReluBiasAdd:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
Relu
IdentityIdentityRelu:activations:0^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*
T0*'
_output_shapes
:џџџџџџџџџ2

Identity"
identityIdentity:output:0*.
_input_shapes
:џџџџџџџџџ::20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:O K
'
_output_shapes
:џџџџџџџџџ
 
_user_specified_nameinputs

з
'__inference_cgan_31_layer_call_fn_71487
input_1
unknown
	unknown_0
	unknown_1
	unknown_2
	unknown_3
	unknown_4
	unknown_5
	unknown_6
identityЂStatefulPartitionedCallС
StatefulPartitionedCallStatefulPartitionedCallinput_1unknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6*
Tin
2	*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:џџџџџџџџџ**
_read_only_resource_inputs

*-
config_proto

CPU

GPU 2J 8 *K
fFRD
B__inference_cgan_31_layer_call_and_return_conditional_losses_713312
StatefulPartitionedCall
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*'
_output_shapes
:џџџџџџџџџ2

Identity"
identityIdentity:output:0*F
_input_shapes5
3:џџџџџџџџџ::::::::22
StatefulPartitionedCallStatefulPartitionedCall:P L
'
_output_shapes
:џџџџџџџџџ
!
_user_specified_name	input_1
+
Ћ
D__inference_generator_layer_call_and_return_conditional_losses_71625

inputs/
+generator_l1_matmul_readvariableop_resource0
,generator_l1_biasadd_readvariableop_resource/
+generator_l2_matmul_readvariableop_resource0
,generator_l2_biasadd_readvariableop_resource/
+generator_l3_matmul_readvariableop_resource0
,generator_l3_biasadd_readvariableop_resource3
/generator_output_matmul_readvariableop_resource4
0generator_output_biasadd_readvariableop_resource
identityЂ#generator_l1/BiasAdd/ReadVariableOpЂ"generator_l1/MatMul/ReadVariableOpЂ#generator_l2/BiasAdd/ReadVariableOpЂ"generator_l2/MatMul/ReadVariableOpЂ#generator_l3/BiasAdd/ReadVariableOpЂ"generator_l3/MatMul/ReadVariableOpЂ'generator_output/BiasAdd/ReadVariableOpЂ&generator_output/MatMul/ReadVariableOpД
"generator_l1/MatMul/ReadVariableOpReadVariableOp+generator_l1_matmul_readvariableop_resource*
_output_shapes

:
*
dtype02$
"generator_l1/MatMul/ReadVariableOp
generator_l1/MatMulMatMulinputs*generator_l1/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ
2
generator_l1/MatMulГ
#generator_l1/BiasAdd/ReadVariableOpReadVariableOp,generator_l1_biasadd_readvariableop_resource*
_output_shapes
:
*
dtype02%
#generator_l1/BiasAdd/ReadVariableOpЕ
generator_l1/BiasAddBiasAddgenerator_l1/MatMul:product:0+generator_l1/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ
2
generator_l1/BiasAdd
generator_l1/ReluRelugenerator_l1/BiasAdd:output:0*
T0*'
_output_shapes
:џџџџџџџџџ
2
generator_l1/ReluД
"generator_l2/MatMul/ReadVariableOpReadVariableOp+generator_l2_matmul_readvariableop_resource*
_output_shapes

:
*
dtype02$
"generator_l2/MatMul/ReadVariableOpГ
generator_l2/MatMulMatMulgenerator_l1/Relu:activations:0*generator_l2/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2
generator_l2/MatMulГ
#generator_l2/BiasAdd/ReadVariableOpReadVariableOp,generator_l2_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02%
#generator_l2/BiasAdd/ReadVariableOpЕ
generator_l2/BiasAddBiasAddgenerator_l2/MatMul:product:0+generator_l2/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2
generator_l2/BiasAdd
generator_l2/ReluRelugenerator_l2/BiasAdd:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
generator_l2/ReluД
"generator_l3/MatMul/ReadVariableOpReadVariableOp+generator_l3_matmul_readvariableop_resource*
_output_shapes

:*
dtype02$
"generator_l3/MatMul/ReadVariableOpГ
generator_l3/MatMulMatMulgenerator_l2/Relu:activations:0*generator_l3/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2
generator_l3/MatMulГ
#generator_l3/BiasAdd/ReadVariableOpReadVariableOp,generator_l3_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02%
#generator_l3/BiasAdd/ReadVariableOpЕ
generator_l3/BiasAddBiasAddgenerator_l3/MatMul:product:0+generator_l3/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2
generator_l3/BiasAdd
generator_l3/ReluRelugenerator_l3/BiasAdd:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
generator_l3/ReluР
&generator_output/MatMul/ReadVariableOpReadVariableOp/generator_output_matmul_readvariableop_resource*
_output_shapes

:*
dtype02(
&generator_output/MatMul/ReadVariableOpП
generator_output/MatMulMatMulgenerator_l3/Relu:activations:0.generator_output/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2
generator_output/MatMulП
'generator_output/BiasAdd/ReadVariableOpReadVariableOp0generator_output_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02)
'generator_output/BiasAdd/ReadVariableOpХ
generator_output/BiasAddBiasAdd!generator_output/MatMul:product:0/generator_output/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2
generator_output/BiasAdd
generator_output/ReluRelu!generator_output/BiasAdd:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
generator_output/ReluЋ
IdentityIdentity#generator_output/Relu:activations:0$^generator_l1/BiasAdd/ReadVariableOp#^generator_l1/MatMul/ReadVariableOp$^generator_l2/BiasAdd/ReadVariableOp#^generator_l2/MatMul/ReadVariableOp$^generator_l3/BiasAdd/ReadVariableOp#^generator_l3/MatMul/ReadVariableOp(^generator_output/BiasAdd/ReadVariableOp'^generator_output/MatMul/ReadVariableOp*
T0*'
_output_shapes
:џџџџџџџџџ2

Identity"
identityIdentity:output:0*F
_input_shapes5
3:џџџџџџџџџ::::::::2J
#generator_l1/BiasAdd/ReadVariableOp#generator_l1/BiasAdd/ReadVariableOp2H
"generator_l1/MatMul/ReadVariableOp"generator_l1/MatMul/ReadVariableOp2J
#generator_l2/BiasAdd/ReadVariableOp#generator_l2/BiasAdd/ReadVariableOp2H
"generator_l2/MatMul/ReadVariableOp"generator_l2/MatMul/ReadVariableOp2J
#generator_l3/BiasAdd/ReadVariableOp#generator_l3/BiasAdd/ReadVariableOp2H
"generator_l3/MatMul/ReadVariableOp"generator_l3/MatMul/ReadVariableOp2R
'generator_output/BiasAdd/ReadVariableOp'generator_output/BiasAdd/ReadVariableOp2P
&generator_output/MatMul/ReadVariableOp&generator_output/MatMul/ReadVariableOp:O K
'
_output_shapes
:џџџџџџџџџ
 
_user_specified_nameinputs

Ы
D__inference_generator_layer_call_and_return_conditional_losses_70933

inputs
generator_l1_70912
generator_l1_70914
generator_l2_70917
generator_l2_70919
generator_l3_70922
generator_l3_70924
generator_output_70927
generator_output_70929
identityЂ$generator_l1/StatefulPartitionedCallЂ$generator_l2/StatefulPartitionedCallЂ$generator_l3/StatefulPartitionedCallЂ(generator_output/StatefulPartitionedCallЅ
$generator_l1/StatefulPartitionedCallStatefulPartitionedCallinputsgenerator_l1_70912generator_l1_70914*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:џџџџџџџџџ
*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8 *P
fKRI
G__inference_generator_l1_layer_call_and_return_conditional_losses_707392&
$generator_l1/StatefulPartitionedCallЬ
$generator_l2/StatefulPartitionedCallStatefulPartitionedCall-generator_l1/StatefulPartitionedCall:output:0generator_l2_70917generator_l2_70919*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:џџџџџџџџџ*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8 *P
fKRI
G__inference_generator_l2_layer_call_and_return_conditional_losses_707662&
$generator_l2/StatefulPartitionedCallЬ
$generator_l3/StatefulPartitionedCallStatefulPartitionedCall-generator_l2/StatefulPartitionedCall:output:0generator_l3_70922generator_l3_70924*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:џџџџџџџџџ*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8 *P
fKRI
G__inference_generator_l3_layer_call_and_return_conditional_losses_707932&
$generator_l3/StatefulPartitionedCallр
(generator_output/StatefulPartitionedCallStatefulPartitionedCall-generator_l3/StatefulPartitionedCall:output:0generator_output_70927generator_output_70929*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:џџџџџџџџџ*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8 *T
fORM
K__inference_generator_output_layer_call_and_return_conditional_losses_708202*
(generator_output/StatefulPartitionedCallЅ
IdentityIdentity1generator_output/StatefulPartitionedCall:output:0%^generator_l1/StatefulPartitionedCall%^generator_l2/StatefulPartitionedCall%^generator_l3/StatefulPartitionedCall)^generator_output/StatefulPartitionedCall*
T0*'
_output_shapes
:џџџџџџџџџ2

Identity"
identityIdentity:output:0*F
_input_shapes5
3:џџџџџџџџџ::::::::2L
$generator_l1/StatefulPartitionedCall$generator_l1/StatefulPartitionedCall2L
$generator_l2/StatefulPartitionedCall$generator_l2/StatefulPartitionedCall2L
$generator_l3/StatefulPartitionedCall$generator_l3/StatefulPartitionedCall2T
(generator_output/StatefulPartitionedCall(generator_output/StatefulPartitionedCall:O K
'
_output_shapes
:џџџџџџџџџ
 
_user_specified_nameinputs
В:
Ъ
B__inference_cgan_31_layer_call_and_return_conditional_losses_71413
input_1A
=discriminator_discriminator_l1_matmul_readvariableop_resourceB
>discriminator_discriminator_l1_biasadd_readvariableop_resourceA
=discriminator_discriminator_l2_matmul_readvariableop_resourceB
>discriminator_discriminator_l2_biasadd_readvariableop_resourceA
=discriminator_discriminator_l3_matmul_readvariableop_resourceB
>discriminator_discriminator_l3_biasadd_readvariableop_resourceE
Adiscriminator_discriminator_output_matmul_readvariableop_resourceF
Bdiscriminator_discriminator_output_biasadd_readvariableop_resource
identityЂ5discriminator/discriminator_l1/BiasAdd/ReadVariableOpЂ4discriminator/discriminator_l1/MatMul/ReadVariableOpЂ5discriminator/discriminator_l2/BiasAdd/ReadVariableOpЂ4discriminator/discriminator_l2/MatMul/ReadVariableOpЂ5discriminator/discriminator_l3/BiasAdd/ReadVariableOpЂ4discriminator/discriminator_l3/MatMul/ReadVariableOpЂ9discriminator/discriminator_output/BiasAdd/ReadVariableOpЂ8discriminator/discriminator_output/MatMul/ReadVariableOpъ
4discriminator/discriminator_l1/MatMul/ReadVariableOpReadVariableOp=discriminator_discriminator_l1_matmul_readvariableop_resource*
_output_shapes

:
*
dtype026
4discriminator/discriminator_l1/MatMul/ReadVariableOpб
%discriminator/discriminator_l1/MatMulMatMulinput_1<discriminator/discriminator_l1/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ
2'
%discriminator/discriminator_l1/MatMulщ
5discriminator/discriminator_l1/BiasAdd/ReadVariableOpReadVariableOp>discriminator_discriminator_l1_biasadd_readvariableop_resource*
_output_shapes
:
*
dtype027
5discriminator/discriminator_l1/BiasAdd/ReadVariableOp§
&discriminator/discriminator_l1/BiasAddBiasAdd/discriminator/discriminator_l1/MatMul:product:0=discriminator/discriminator_l1/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ
2(
&discriminator/discriminator_l1/BiasAddЕ
#discriminator/discriminator_l1/ReluRelu/discriminator/discriminator_l1/BiasAdd:output:0*
T0*'
_output_shapes
:џџџџџџџџџ
2%
#discriminator/discriminator_l1/Reluъ
4discriminator/discriminator_l2/MatMul/ReadVariableOpReadVariableOp=discriminator_discriminator_l2_matmul_readvariableop_resource*
_output_shapes

:
*
dtype026
4discriminator/discriminator_l2/MatMul/ReadVariableOpћ
%discriminator/discriminator_l2/MatMulMatMul1discriminator/discriminator_l1/Relu:activations:0<discriminator/discriminator_l2/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2'
%discriminator/discriminator_l2/MatMulщ
5discriminator/discriminator_l2/BiasAdd/ReadVariableOpReadVariableOp>discriminator_discriminator_l2_biasadd_readvariableop_resource*
_output_shapes
:*
dtype027
5discriminator/discriminator_l2/BiasAdd/ReadVariableOp§
&discriminator/discriminator_l2/BiasAddBiasAdd/discriminator/discriminator_l2/MatMul:product:0=discriminator/discriminator_l2/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2(
&discriminator/discriminator_l2/BiasAddЕ
#discriminator/discriminator_l2/ReluRelu/discriminator/discriminator_l2/BiasAdd:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2%
#discriminator/discriminator_l2/Reluъ
4discriminator/discriminator_l3/MatMul/ReadVariableOpReadVariableOp=discriminator_discriminator_l3_matmul_readvariableop_resource*
_output_shapes

:*
dtype026
4discriminator/discriminator_l3/MatMul/ReadVariableOpћ
%discriminator/discriminator_l3/MatMulMatMul1discriminator/discriminator_l2/Relu:activations:0<discriminator/discriminator_l3/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2'
%discriminator/discriminator_l3/MatMulщ
5discriminator/discriminator_l3/BiasAdd/ReadVariableOpReadVariableOp>discriminator_discriminator_l3_biasadd_readvariableop_resource*
_output_shapes
:*
dtype027
5discriminator/discriminator_l3/BiasAdd/ReadVariableOp§
&discriminator/discriminator_l3/BiasAddBiasAdd/discriminator/discriminator_l3/MatMul:product:0=discriminator/discriminator_l3/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2(
&discriminator/discriminator_l3/BiasAddЕ
#discriminator/discriminator_l3/ReluRelu/discriminator/discriminator_l3/BiasAdd:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2%
#discriminator/discriminator_l3/Reluі
8discriminator/discriminator_output/MatMul/ReadVariableOpReadVariableOpAdiscriminator_discriminator_output_matmul_readvariableop_resource*
_output_shapes

:*
dtype02:
8discriminator/discriminator_output/MatMul/ReadVariableOp
)discriminator/discriminator_output/MatMulMatMul1discriminator/discriminator_l3/Relu:activations:0@discriminator/discriminator_output/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2+
)discriminator/discriminator_output/MatMulѕ
9discriminator/discriminator_output/BiasAdd/ReadVariableOpReadVariableOpBdiscriminator_discriminator_output_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02;
9discriminator/discriminator_output/BiasAdd/ReadVariableOp
*discriminator/discriminator_output/BiasAddBiasAdd3discriminator/discriminator_output/MatMul:product:0Adiscriminator/discriminator_output/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2,
*discriminator/discriminator_output/BiasAddЪ
*discriminator/discriminator_output/SigmoidSigmoid3discriminator/discriminator_output/BiasAdd:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2,
*discriminator/discriminator_output/SigmoidЦ
IdentityIdentity.discriminator/discriminator_output/Sigmoid:y:06^discriminator/discriminator_l1/BiasAdd/ReadVariableOp5^discriminator/discriminator_l1/MatMul/ReadVariableOp6^discriminator/discriminator_l2/BiasAdd/ReadVariableOp5^discriminator/discriminator_l2/MatMul/ReadVariableOp6^discriminator/discriminator_l3/BiasAdd/ReadVariableOp5^discriminator/discriminator_l3/MatMul/ReadVariableOp:^discriminator/discriminator_output/BiasAdd/ReadVariableOp9^discriminator/discriminator_output/MatMul/ReadVariableOp*
T0*'
_output_shapes
:џџџџџџџџџ2

Identity"
identityIdentity:output:0*F
_input_shapes5
3:џџџџџџџџџ::::::::2n
5discriminator/discriminator_l1/BiasAdd/ReadVariableOp5discriminator/discriminator_l1/BiasAdd/ReadVariableOp2l
4discriminator/discriminator_l1/MatMul/ReadVariableOp4discriminator/discriminator_l1/MatMul/ReadVariableOp2n
5discriminator/discriminator_l2/BiasAdd/ReadVariableOp5discriminator/discriminator_l2/BiasAdd/ReadVariableOp2l
4discriminator/discriminator_l2/MatMul/ReadVariableOp4discriminator/discriminator_l2/MatMul/ReadVariableOp2n
5discriminator/discriminator_l3/BiasAdd/ReadVariableOp5discriminator/discriminator_l3/BiasAdd/ReadVariableOp2l
4discriminator/discriminator_l3/MatMul/ReadVariableOp4discriminator/discriminator_l3/MatMul/ReadVariableOp2v
9discriminator/discriminator_output/BiasAdd/ReadVariableOp9discriminator/discriminator_output/BiasAdd/ReadVariableOp2t
8discriminator/discriminator_output/MatMul/ReadVariableOp8discriminator/discriminator_output/MatMul/ReadVariableOp:P L
'
_output_shapes
:џџџџџџџџџ
!
_user_specified_name	input_1

и
)__inference_generator_layer_call_fn_71678

inputs
unknown
	unknown_0
	unknown_1
	unknown_2
	unknown_3
	unknown_4
	unknown_5
	unknown_6
identityЂStatefulPartitionedCallТ
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6*
Tin
2	*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:џџџџџџџџџ**
_read_only_resource_inputs

*-
config_proto

CPU

GPU 2J 8 *M
fHRF
D__inference_generator_layer_call_and_return_conditional_losses_708882
StatefulPartitionedCall
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*'
_output_shapes
:џџџџџџџџџ2

Identity"
identityIdentity:output:0*F
_input_shapes5
3:џџџџџџџџџ::::::::22
StatefulPartitionedCallStatefulPartitionedCall:O K
'
_output_shapes
:џџџџџџџџџ
 
_user_specified_nameinputs
Р

б
B__inference_cgan_31_layer_call_and_return_conditional_losses_71289

inputs
discriminator_71271
discriminator_71273
discriminator_71275
discriminator_71277
discriminator_71279
discriminator_71281
discriminator_71283
discriminator_71285
identityЂ%discriminator/StatefulPartitionedCallД
%discriminator/StatefulPartitionedCallStatefulPartitionedCallinputsdiscriminator_71271discriminator_71273discriminator_71275discriminator_71277discriminator_71279discriminator_71281discriminator_71283discriminator_71285*
Tin
2	*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:џџџџџџџџџ**
_read_only_resource_inputs

*-
config_proto

CPU

GPU 2J 8 *Q
fLRJ
H__inference_discriminator_layer_call_and_return_conditional_losses_711162'
%discriminator/StatefulPartitionedCallЊ
IdentityIdentity.discriminator/StatefulPartitionedCall:output:0&^discriminator/StatefulPartitionedCall*
T0*'
_output_shapes
:џџџџџџџџџ2

Identity"
identityIdentity:output:0*F
_input_shapes5
3:џџџџџџџџџ::::::::2N
%discriminator/StatefulPartitionedCall%discriminator/StatefulPartitionedCall:O K
'
_output_shapes
:џџџџџџџџџ
 
_user_specified_nameinputs
В
д
D__inference_generator_layer_call_and_return_conditional_losses_70861
generator_input
generator_l1_70840
generator_l1_70842
generator_l2_70845
generator_l2_70847
generator_l3_70850
generator_l3_70852
generator_output_70855
generator_output_70857
identityЂ$generator_l1/StatefulPartitionedCallЂ$generator_l2/StatefulPartitionedCallЂ$generator_l3/StatefulPartitionedCallЂ(generator_output/StatefulPartitionedCallЎ
$generator_l1/StatefulPartitionedCallStatefulPartitionedCallgenerator_inputgenerator_l1_70840generator_l1_70842*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:џџџџџџџџџ
*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8 *P
fKRI
G__inference_generator_l1_layer_call_and_return_conditional_losses_707392&
$generator_l1/StatefulPartitionedCallЬ
$generator_l2/StatefulPartitionedCallStatefulPartitionedCall-generator_l1/StatefulPartitionedCall:output:0generator_l2_70845generator_l2_70847*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:џџџџџџџџџ*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8 *P
fKRI
G__inference_generator_l2_layer_call_and_return_conditional_losses_707662&
$generator_l2/StatefulPartitionedCallЬ
$generator_l3/StatefulPartitionedCallStatefulPartitionedCall-generator_l2/StatefulPartitionedCall:output:0generator_l3_70850generator_l3_70852*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:џџџџџџџџџ*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8 *P
fKRI
G__inference_generator_l3_layer_call_and_return_conditional_losses_707932&
$generator_l3/StatefulPartitionedCallр
(generator_output/StatefulPartitionedCallStatefulPartitionedCall-generator_l3/StatefulPartitionedCall:output:0generator_output_70855generator_output_70857*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:џџџџџџџџџ*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8 *T
fORM
K__inference_generator_output_layer_call_and_return_conditional_losses_708202*
(generator_output/StatefulPartitionedCallЅ
IdentityIdentity1generator_output/StatefulPartitionedCall:output:0%^generator_l1/StatefulPartitionedCall%^generator_l2/StatefulPartitionedCall%^generator_l3/StatefulPartitionedCall)^generator_output/StatefulPartitionedCall*
T0*'
_output_shapes
:џџџџџџџџџ2

Identity"
identityIdentity:output:0*F
_input_shapes5
3:џџџџџџџџџ::::::::2L
$generator_l1/StatefulPartitionedCall$generator_l1/StatefulPartitionedCall2L
$generator_l2/StatefulPartitionedCall$generator_l2/StatefulPartitionedCall2L
$generator_l3/StatefulPartitionedCall$generator_l3/StatefulPartitionedCall2T
(generator_output/StatefulPartitionedCall(generator_output/StatefulPartitionedCall:X T
'
_output_shapes
:џџџџџџџџџ
)
_user_specified_namegenerator_input
В:
Ъ
B__inference_cgan_31_layer_call_and_return_conditional_losses_71445
input_1A
=discriminator_discriminator_l1_matmul_readvariableop_resourceB
>discriminator_discriminator_l1_biasadd_readvariableop_resourceA
=discriminator_discriminator_l2_matmul_readvariableop_resourceB
>discriminator_discriminator_l2_biasadd_readvariableop_resourceA
=discriminator_discriminator_l3_matmul_readvariableop_resourceB
>discriminator_discriminator_l3_biasadd_readvariableop_resourceE
Adiscriminator_discriminator_output_matmul_readvariableop_resourceF
Bdiscriminator_discriminator_output_biasadd_readvariableop_resource
identityЂ5discriminator/discriminator_l1/BiasAdd/ReadVariableOpЂ4discriminator/discriminator_l1/MatMul/ReadVariableOpЂ5discriminator/discriminator_l2/BiasAdd/ReadVariableOpЂ4discriminator/discriminator_l2/MatMul/ReadVariableOpЂ5discriminator/discriminator_l3/BiasAdd/ReadVariableOpЂ4discriminator/discriminator_l3/MatMul/ReadVariableOpЂ9discriminator/discriminator_output/BiasAdd/ReadVariableOpЂ8discriminator/discriminator_output/MatMul/ReadVariableOpъ
4discriminator/discriminator_l1/MatMul/ReadVariableOpReadVariableOp=discriminator_discriminator_l1_matmul_readvariableop_resource*
_output_shapes

:
*
dtype026
4discriminator/discriminator_l1/MatMul/ReadVariableOpб
%discriminator/discriminator_l1/MatMulMatMulinput_1<discriminator/discriminator_l1/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ
2'
%discriminator/discriminator_l1/MatMulщ
5discriminator/discriminator_l1/BiasAdd/ReadVariableOpReadVariableOp>discriminator_discriminator_l1_biasadd_readvariableop_resource*
_output_shapes
:
*
dtype027
5discriminator/discriminator_l1/BiasAdd/ReadVariableOp§
&discriminator/discriminator_l1/BiasAddBiasAdd/discriminator/discriminator_l1/MatMul:product:0=discriminator/discriminator_l1/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ
2(
&discriminator/discriminator_l1/BiasAddЕ
#discriminator/discriminator_l1/ReluRelu/discriminator/discriminator_l1/BiasAdd:output:0*
T0*'
_output_shapes
:џџџџџџџџџ
2%
#discriminator/discriminator_l1/Reluъ
4discriminator/discriminator_l2/MatMul/ReadVariableOpReadVariableOp=discriminator_discriminator_l2_matmul_readvariableop_resource*
_output_shapes

:
*
dtype026
4discriminator/discriminator_l2/MatMul/ReadVariableOpћ
%discriminator/discriminator_l2/MatMulMatMul1discriminator/discriminator_l1/Relu:activations:0<discriminator/discriminator_l2/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2'
%discriminator/discriminator_l2/MatMulщ
5discriminator/discriminator_l2/BiasAdd/ReadVariableOpReadVariableOp>discriminator_discriminator_l2_biasadd_readvariableop_resource*
_output_shapes
:*
dtype027
5discriminator/discriminator_l2/BiasAdd/ReadVariableOp§
&discriminator/discriminator_l2/BiasAddBiasAdd/discriminator/discriminator_l2/MatMul:product:0=discriminator/discriminator_l2/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2(
&discriminator/discriminator_l2/BiasAddЕ
#discriminator/discriminator_l2/ReluRelu/discriminator/discriminator_l2/BiasAdd:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2%
#discriminator/discriminator_l2/Reluъ
4discriminator/discriminator_l3/MatMul/ReadVariableOpReadVariableOp=discriminator_discriminator_l3_matmul_readvariableop_resource*
_output_shapes

:*
dtype026
4discriminator/discriminator_l3/MatMul/ReadVariableOpћ
%discriminator/discriminator_l3/MatMulMatMul1discriminator/discriminator_l2/Relu:activations:0<discriminator/discriminator_l3/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2'
%discriminator/discriminator_l3/MatMulщ
5discriminator/discriminator_l3/BiasAdd/ReadVariableOpReadVariableOp>discriminator_discriminator_l3_biasadd_readvariableop_resource*
_output_shapes
:*
dtype027
5discriminator/discriminator_l3/BiasAdd/ReadVariableOp§
&discriminator/discriminator_l3/BiasAddBiasAdd/discriminator/discriminator_l3/MatMul:product:0=discriminator/discriminator_l3/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2(
&discriminator/discriminator_l3/BiasAddЕ
#discriminator/discriminator_l3/ReluRelu/discriminator/discriminator_l3/BiasAdd:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2%
#discriminator/discriminator_l3/Reluі
8discriminator/discriminator_output/MatMul/ReadVariableOpReadVariableOpAdiscriminator_discriminator_output_matmul_readvariableop_resource*
_output_shapes

:*
dtype02:
8discriminator/discriminator_output/MatMul/ReadVariableOp
)discriminator/discriminator_output/MatMulMatMul1discriminator/discriminator_l3/Relu:activations:0@discriminator/discriminator_output/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2+
)discriminator/discriminator_output/MatMulѕ
9discriminator/discriminator_output/BiasAdd/ReadVariableOpReadVariableOpBdiscriminator_discriminator_output_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02;
9discriminator/discriminator_output/BiasAdd/ReadVariableOp
*discriminator/discriminator_output/BiasAddBiasAdd3discriminator/discriminator_output/MatMul:product:0Adiscriminator/discriminator_output/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2,
*discriminator/discriminator_output/BiasAddЪ
*discriminator/discriminator_output/SigmoidSigmoid3discriminator/discriminator_output/BiasAdd:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2,
*discriminator/discriminator_output/SigmoidЦ
IdentityIdentity.discriminator/discriminator_output/Sigmoid:y:06^discriminator/discriminator_l1/BiasAdd/ReadVariableOp5^discriminator/discriminator_l1/MatMul/ReadVariableOp6^discriminator/discriminator_l2/BiasAdd/ReadVariableOp5^discriminator/discriminator_l2/MatMul/ReadVariableOp6^discriminator/discriminator_l3/BiasAdd/ReadVariableOp5^discriminator/discriminator_l3/MatMul/ReadVariableOp:^discriminator/discriminator_output/BiasAdd/ReadVariableOp9^discriminator/discriminator_output/MatMul/ReadVariableOp*
T0*'
_output_shapes
:џџџџџџџџџ2

Identity"
identityIdentity:output:0*F
_input_shapes5
3:џџџџџџџџџ::::::::2n
5discriminator/discriminator_l1/BiasAdd/ReadVariableOp5discriminator/discriminator_l1/BiasAdd/ReadVariableOp2l
4discriminator/discriminator_l1/MatMul/ReadVariableOp4discriminator/discriminator_l1/MatMul/ReadVariableOp2n
5discriminator/discriminator_l2/BiasAdd/ReadVariableOp5discriminator/discriminator_l2/BiasAdd/ReadVariableOp2l
4discriminator/discriminator_l2/MatMul/ReadVariableOp4discriminator/discriminator_l2/MatMul/ReadVariableOp2n
5discriminator/discriminator_l3/BiasAdd/ReadVariableOp5discriminator/discriminator_l3/BiasAdd/ReadVariableOp2l
4discriminator/discriminator_l3/MatMul/ReadVariableOp4discriminator/discriminator_l3/MatMul/ReadVariableOp2v
9discriminator/discriminator_output/BiasAdd/ReadVariableOp9discriminator/discriminator_output/BiasAdd/ReadVariableOp2t
8discriminator/discriminator_output/MatMul/ReadVariableOp8discriminator/discriminator_output/MatMul/ReadVariableOp:P L
'
_output_shapes
:џџџџџџџџџ
!
_user_specified_name	input_1
ѕ	
ф
K__inference_discriminator_l2_layer_call_and_return_conditional_losses_70994

inputs"
matmul_readvariableop_resource#
biasadd_readvariableop_resource
identityЂBiasAdd/ReadVariableOpЂMatMul/ReadVariableOp
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:
*
dtype02
MatMul/ReadVariableOps
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2
MatMul
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:*
dtype02
BiasAdd/ReadVariableOp
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2	
BiasAddX
ReluReluBiasAdd:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
Relu
IdentityIdentityRelu:activations:0^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*
T0*'
_output_shapes
:џџџџџџџџџ2

Identity"
identityIdentity:output:0*.
_input_shapes
:џџџџџџџџџ
::20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:O K
'
_output_shapes
:џџџџџџџџџ

 
_user_specified_nameinputs
ћ	
ш
O__inference_discriminator_output_layer_call_and_return_conditional_losses_71048

inputs"
matmul_readvariableop_resource#
biasadd_readvariableop_resource
identityЂBiasAdd/ReadVariableOpЂMatMul/ReadVariableOp
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:*
dtype02
MatMul/ReadVariableOps
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2
MatMul
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:*
dtype02
BiasAdd/ReadVariableOp
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2	
BiasAdda
SigmoidSigmoidBiasAdd:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2	
Sigmoid
IdentityIdentitySigmoid:y:0^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*
T0*'
_output_shapes
:џџџџџџџџџ2

Identity"
identityIdentity:output:0*.
_input_shapes
:џџџџџџџџџ::20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:O K
'
_output_shapes
:џџџџџџџџџ
 
_user_specified_nameinputs
ѓ

4__inference_discriminator_output_layer_call_fn_71965

inputs
unknown
	unknown_0
identityЂStatefulPartitionedCallџ
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:џџџџџџџџџ*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8 *X
fSRQ
O__inference_discriminator_output_layer_call_and_return_conditional_losses_710482
StatefulPartitionedCall
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*'
_output_shapes
:џџџџџџџџџ2

Identity"
identityIdentity:output:0*.
_input_shapes
:џџџџџџџџџ::22
StatefulPartitionedCallStatefulPartitionedCall:O K
'
_output_shapes
:џџџџџџџџџ
 
_user_specified_nameinputs
Э.
я
H__inference_discriminator_layer_call_and_return_conditional_losses_71731

inputs3
/discriminator_l1_matmul_readvariableop_resource4
0discriminator_l1_biasadd_readvariableop_resource3
/discriminator_l2_matmul_readvariableop_resource4
0discriminator_l2_biasadd_readvariableop_resource3
/discriminator_l3_matmul_readvariableop_resource4
0discriminator_l3_biasadd_readvariableop_resource7
3discriminator_output_matmul_readvariableop_resource8
4discriminator_output_biasadd_readvariableop_resource
identityЂ'discriminator_l1/BiasAdd/ReadVariableOpЂ&discriminator_l1/MatMul/ReadVariableOpЂ'discriminator_l2/BiasAdd/ReadVariableOpЂ&discriminator_l2/MatMul/ReadVariableOpЂ'discriminator_l3/BiasAdd/ReadVariableOpЂ&discriminator_l3/MatMul/ReadVariableOpЂ+discriminator_output/BiasAdd/ReadVariableOpЂ*discriminator_output/MatMul/ReadVariableOpР
&discriminator_l1/MatMul/ReadVariableOpReadVariableOp/discriminator_l1_matmul_readvariableop_resource*
_output_shapes

:
*
dtype02(
&discriminator_l1/MatMul/ReadVariableOpІ
discriminator_l1/MatMulMatMulinputs.discriminator_l1/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ
2
discriminator_l1/MatMulП
'discriminator_l1/BiasAdd/ReadVariableOpReadVariableOp0discriminator_l1_biasadd_readvariableop_resource*
_output_shapes
:
*
dtype02)
'discriminator_l1/BiasAdd/ReadVariableOpХ
discriminator_l1/BiasAddBiasAdd!discriminator_l1/MatMul:product:0/discriminator_l1/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ
2
discriminator_l1/BiasAdd
discriminator_l1/ReluRelu!discriminator_l1/BiasAdd:output:0*
T0*'
_output_shapes
:џџџџџџџџџ
2
discriminator_l1/ReluР
&discriminator_l2/MatMul/ReadVariableOpReadVariableOp/discriminator_l2_matmul_readvariableop_resource*
_output_shapes

:
*
dtype02(
&discriminator_l2/MatMul/ReadVariableOpУ
discriminator_l2/MatMulMatMul#discriminator_l1/Relu:activations:0.discriminator_l2/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2
discriminator_l2/MatMulП
'discriminator_l2/BiasAdd/ReadVariableOpReadVariableOp0discriminator_l2_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02)
'discriminator_l2/BiasAdd/ReadVariableOpХ
discriminator_l2/BiasAddBiasAdd!discriminator_l2/MatMul:product:0/discriminator_l2/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2
discriminator_l2/BiasAdd
discriminator_l2/ReluRelu!discriminator_l2/BiasAdd:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
discriminator_l2/ReluР
&discriminator_l3/MatMul/ReadVariableOpReadVariableOp/discriminator_l3_matmul_readvariableop_resource*
_output_shapes

:*
dtype02(
&discriminator_l3/MatMul/ReadVariableOpУ
discriminator_l3/MatMulMatMul#discriminator_l2/Relu:activations:0.discriminator_l3/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2
discriminator_l3/MatMulП
'discriminator_l3/BiasAdd/ReadVariableOpReadVariableOp0discriminator_l3_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02)
'discriminator_l3/BiasAdd/ReadVariableOpХ
discriminator_l3/BiasAddBiasAdd!discriminator_l3/MatMul:product:0/discriminator_l3/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2
discriminator_l3/BiasAdd
discriminator_l3/ReluRelu!discriminator_l3/BiasAdd:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
discriminator_l3/ReluЬ
*discriminator_output/MatMul/ReadVariableOpReadVariableOp3discriminator_output_matmul_readvariableop_resource*
_output_shapes

:*
dtype02,
*discriminator_output/MatMul/ReadVariableOpЯ
discriminator_output/MatMulMatMul#discriminator_l3/Relu:activations:02discriminator_output/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2
discriminator_output/MatMulЫ
+discriminator_output/BiasAdd/ReadVariableOpReadVariableOp4discriminator_output_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02-
+discriminator_output/BiasAdd/ReadVariableOpе
discriminator_output/BiasAddBiasAdd%discriminator_output/MatMul:product:03discriminator_output/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2
discriminator_output/BiasAdd 
discriminator_output/SigmoidSigmoid%discriminator_output/BiasAdd:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
discriminator_output/SigmoidШ
IdentityIdentity discriminator_output/Sigmoid:y:0(^discriminator_l1/BiasAdd/ReadVariableOp'^discriminator_l1/MatMul/ReadVariableOp(^discriminator_l2/BiasAdd/ReadVariableOp'^discriminator_l2/MatMul/ReadVariableOp(^discriminator_l3/BiasAdd/ReadVariableOp'^discriminator_l3/MatMul/ReadVariableOp,^discriminator_output/BiasAdd/ReadVariableOp+^discriminator_output/MatMul/ReadVariableOp*
T0*'
_output_shapes
:џџџџџџџџџ2

Identity"
identityIdentity:output:0*F
_input_shapes5
3:џџџџџџџџџ::::::::2R
'discriminator_l1/BiasAdd/ReadVariableOp'discriminator_l1/BiasAdd/ReadVariableOp2P
&discriminator_l1/MatMul/ReadVariableOp&discriminator_l1/MatMul/ReadVariableOp2R
'discriminator_l2/BiasAdd/ReadVariableOp'discriminator_l2/BiasAdd/ReadVariableOp2P
&discriminator_l2/MatMul/ReadVariableOp&discriminator_l2/MatMul/ReadVariableOp2R
'discriminator_l3/BiasAdd/ReadVariableOp'discriminator_l3/BiasAdd/ReadVariableOp2P
&discriminator_l3/MatMul/ReadVariableOp&discriminator_l3/MatMul/ReadVariableOp2Z
+discriminator_output/BiasAdd/ReadVariableOp+discriminator_output/BiasAdd/ReadVariableOp2X
*discriminator_output/MatMul/ReadVariableOp*discriminator_output/MatMul/ReadVariableOp:O K
'
_output_shapes
:џџџџџџџџџ
 
_user_specified_nameinputs
Э.
я
H__inference_discriminator_layer_call_and_return_conditional_losses_71763

inputs3
/discriminator_l1_matmul_readvariableop_resource4
0discriminator_l1_biasadd_readvariableop_resource3
/discriminator_l2_matmul_readvariableop_resource4
0discriminator_l2_biasadd_readvariableop_resource3
/discriminator_l3_matmul_readvariableop_resource4
0discriminator_l3_biasadd_readvariableop_resource7
3discriminator_output_matmul_readvariableop_resource8
4discriminator_output_biasadd_readvariableop_resource
identityЂ'discriminator_l1/BiasAdd/ReadVariableOpЂ&discriminator_l1/MatMul/ReadVariableOpЂ'discriminator_l2/BiasAdd/ReadVariableOpЂ&discriminator_l2/MatMul/ReadVariableOpЂ'discriminator_l3/BiasAdd/ReadVariableOpЂ&discriminator_l3/MatMul/ReadVariableOpЂ+discriminator_output/BiasAdd/ReadVariableOpЂ*discriminator_output/MatMul/ReadVariableOpР
&discriminator_l1/MatMul/ReadVariableOpReadVariableOp/discriminator_l1_matmul_readvariableop_resource*
_output_shapes

:
*
dtype02(
&discriminator_l1/MatMul/ReadVariableOpІ
discriminator_l1/MatMulMatMulinputs.discriminator_l1/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ
2
discriminator_l1/MatMulП
'discriminator_l1/BiasAdd/ReadVariableOpReadVariableOp0discriminator_l1_biasadd_readvariableop_resource*
_output_shapes
:
*
dtype02)
'discriminator_l1/BiasAdd/ReadVariableOpХ
discriminator_l1/BiasAddBiasAdd!discriminator_l1/MatMul:product:0/discriminator_l1/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ
2
discriminator_l1/BiasAdd
discriminator_l1/ReluRelu!discriminator_l1/BiasAdd:output:0*
T0*'
_output_shapes
:џџџџџџџџџ
2
discriminator_l1/ReluР
&discriminator_l2/MatMul/ReadVariableOpReadVariableOp/discriminator_l2_matmul_readvariableop_resource*
_output_shapes

:
*
dtype02(
&discriminator_l2/MatMul/ReadVariableOpУ
discriminator_l2/MatMulMatMul#discriminator_l1/Relu:activations:0.discriminator_l2/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2
discriminator_l2/MatMulП
'discriminator_l2/BiasAdd/ReadVariableOpReadVariableOp0discriminator_l2_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02)
'discriminator_l2/BiasAdd/ReadVariableOpХ
discriminator_l2/BiasAddBiasAdd!discriminator_l2/MatMul:product:0/discriminator_l2/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2
discriminator_l2/BiasAdd
discriminator_l2/ReluRelu!discriminator_l2/BiasAdd:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
discriminator_l2/ReluР
&discriminator_l3/MatMul/ReadVariableOpReadVariableOp/discriminator_l3_matmul_readvariableop_resource*
_output_shapes

:*
dtype02(
&discriminator_l3/MatMul/ReadVariableOpУ
discriminator_l3/MatMulMatMul#discriminator_l2/Relu:activations:0.discriminator_l3/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2
discriminator_l3/MatMulП
'discriminator_l3/BiasAdd/ReadVariableOpReadVariableOp0discriminator_l3_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02)
'discriminator_l3/BiasAdd/ReadVariableOpХ
discriminator_l3/BiasAddBiasAdd!discriminator_l3/MatMul:product:0/discriminator_l3/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2
discriminator_l3/BiasAdd
discriminator_l3/ReluRelu!discriminator_l3/BiasAdd:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
discriminator_l3/ReluЬ
*discriminator_output/MatMul/ReadVariableOpReadVariableOp3discriminator_output_matmul_readvariableop_resource*
_output_shapes

:*
dtype02,
*discriminator_output/MatMul/ReadVariableOpЯ
discriminator_output/MatMulMatMul#discriminator_l3/Relu:activations:02discriminator_output/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2
discriminator_output/MatMulЫ
+discriminator_output/BiasAdd/ReadVariableOpReadVariableOp4discriminator_output_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02-
+discriminator_output/BiasAdd/ReadVariableOpе
discriminator_output/BiasAddBiasAdd%discriminator_output/MatMul:product:03discriminator_output/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2
discriminator_output/BiasAdd 
discriminator_output/SigmoidSigmoid%discriminator_output/BiasAdd:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
discriminator_output/SigmoidШ
IdentityIdentity discriminator_output/Sigmoid:y:0(^discriminator_l1/BiasAdd/ReadVariableOp'^discriminator_l1/MatMul/ReadVariableOp(^discriminator_l2/BiasAdd/ReadVariableOp'^discriminator_l2/MatMul/ReadVariableOp(^discriminator_l3/BiasAdd/ReadVariableOp'^discriminator_l3/MatMul/ReadVariableOp,^discriminator_output/BiasAdd/ReadVariableOp+^discriminator_output/MatMul/ReadVariableOp*
T0*'
_output_shapes
:џџџџџџџџџ2

Identity"
identityIdentity:output:0*F
_input_shapes5
3:џџџџџџџџџ::::::::2R
'discriminator_l1/BiasAdd/ReadVariableOp'discriminator_l1/BiasAdd/ReadVariableOp2P
&discriminator_l1/MatMul/ReadVariableOp&discriminator_l1/MatMul/ReadVariableOp2R
'discriminator_l2/BiasAdd/ReadVariableOp'discriminator_l2/BiasAdd/ReadVariableOp2P
&discriminator_l2/MatMul/ReadVariableOp&discriminator_l2/MatMul/ReadVariableOp2R
'discriminator_l3/BiasAdd/ReadVariableOp'discriminator_l3/BiasAdd/ReadVariableOp2P
&discriminator_l3/MatMul/ReadVariableOp&discriminator_l3/MatMul/ReadVariableOp2Z
+discriminator_output/BiasAdd/ReadVariableOp+discriminator_output/BiasAdd/ReadVariableOp2X
*discriminator_output/MatMul/ReadVariableOp*discriminator_output/MatMul/ReadVariableOp:O K
'
_output_shapes
:џџџџџџџџџ
 
_user_specified_nameinputs
ё	
р
G__inference_generator_l2_layer_call_and_return_conditional_losses_71836

inputs"
matmul_readvariableop_resource#
biasadd_readvariableop_resource
identityЂBiasAdd/ReadVariableOpЂMatMul/ReadVariableOp
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:
*
dtype02
MatMul/ReadVariableOps
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2
MatMul
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:*
dtype02
BiasAdd/ReadVariableOp
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2	
BiasAddX
ReluReluBiasAdd:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
Relu
IdentityIdentityRelu:activations:0^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*
T0*'
_output_shapes
:џџџџџџџџџ2

Identity"
identityIdentity:output:0*.
_input_shapes
:џџџџџџџџџ
::20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:O K
'
_output_shapes
:џџџџџџџџџ

 
_user_specified_nameinputs
л
џ
H__inference_discriminator_layer_call_and_return_conditional_losses_71161

inputs
discriminator_l1_71140
discriminator_l1_71142
discriminator_l2_71145
discriminator_l2_71147
discriminator_l3_71150
discriminator_l3_71152
discriminator_output_71155
discriminator_output_71157
identityЂ(discriminator_l1/StatefulPartitionedCallЂ(discriminator_l2/StatefulPartitionedCallЂ(discriminator_l3/StatefulPartitionedCallЂ,discriminator_output/StatefulPartitionedCallЙ
(discriminator_l1/StatefulPartitionedCallStatefulPartitionedCallinputsdiscriminator_l1_71140discriminator_l1_71142*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:џџџџџџџџџ
*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8 *T
fORM
K__inference_discriminator_l1_layer_call_and_return_conditional_losses_709672*
(discriminator_l1/StatefulPartitionedCallф
(discriminator_l2/StatefulPartitionedCallStatefulPartitionedCall1discriminator_l1/StatefulPartitionedCall:output:0discriminator_l2_71145discriminator_l2_71147*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:џџџџџџџџџ*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8 *T
fORM
K__inference_discriminator_l2_layer_call_and_return_conditional_losses_709942*
(discriminator_l2/StatefulPartitionedCallф
(discriminator_l3/StatefulPartitionedCallStatefulPartitionedCall1discriminator_l2/StatefulPartitionedCall:output:0discriminator_l3_71150discriminator_l3_71152*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:џџџџџџџџџ*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8 *T
fORM
K__inference_discriminator_l3_layer_call_and_return_conditional_losses_710212*
(discriminator_l3/StatefulPartitionedCallј
,discriminator_output/StatefulPartitionedCallStatefulPartitionedCall1discriminator_l3/StatefulPartitionedCall:output:0discriminator_output_71155discriminator_output_71157*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:џџџџџџџџџ*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8 *X
fSRQ
O__inference_discriminator_output_layer_call_and_return_conditional_losses_710482.
,discriminator_output/StatefulPartitionedCallЙ
IdentityIdentity5discriminator_output/StatefulPartitionedCall:output:0)^discriminator_l1/StatefulPartitionedCall)^discriminator_l2/StatefulPartitionedCall)^discriminator_l3/StatefulPartitionedCall-^discriminator_output/StatefulPartitionedCall*
T0*'
_output_shapes
:џџџџџџџџџ2

Identity"
identityIdentity:output:0*F
_input_shapes5
3:џџџџџџџџџ::::::::2T
(discriminator_l1/StatefulPartitionedCall(discriminator_l1/StatefulPartitionedCall2T
(discriminator_l2/StatefulPartitionedCall(discriminator_l2/StatefulPartitionedCall2T
(discriminator_l3/StatefulPartitionedCall(discriminator_l3/StatefulPartitionedCall2\
,discriminator_output/StatefulPartitionedCall,discriminator_output/StatefulPartitionedCall:O K
'
_output_shapes
:џџџџџџџџџ
 
_user_specified_nameinputs
ѕ	
ф
K__inference_discriminator_l2_layer_call_and_return_conditional_losses_71916

inputs"
matmul_readvariableop_resource#
biasadd_readvariableop_resource
identityЂBiasAdd/ReadVariableOpЂMatMul/ReadVariableOp
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:
*
dtype02
MatMul/ReadVariableOps
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2
MatMul
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:*
dtype02
BiasAdd/ReadVariableOp
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2	
BiasAddX
ReluReluBiasAdd:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
Relu
IdentityIdentityRelu:activations:0^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*
T0*'
_output_shapes
:џџџџџџџџџ2

Identity"
identityIdentity:output:0*.
_input_shapes
:џџџџџџџџџ
::20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:O K
'
_output_shapes
:џџџџџџџџџ

 
_user_specified_nameinputs
Џ:
Щ
B__inference_cgan_31_layer_call_and_return_conditional_losses_71519

inputsA
=discriminator_discriminator_l1_matmul_readvariableop_resourceB
>discriminator_discriminator_l1_biasadd_readvariableop_resourceA
=discriminator_discriminator_l2_matmul_readvariableop_resourceB
>discriminator_discriminator_l2_biasadd_readvariableop_resourceA
=discriminator_discriminator_l3_matmul_readvariableop_resourceB
>discriminator_discriminator_l3_biasadd_readvariableop_resourceE
Adiscriminator_discriminator_output_matmul_readvariableop_resourceF
Bdiscriminator_discriminator_output_biasadd_readvariableop_resource
identityЂ5discriminator/discriminator_l1/BiasAdd/ReadVariableOpЂ4discriminator/discriminator_l1/MatMul/ReadVariableOpЂ5discriminator/discriminator_l2/BiasAdd/ReadVariableOpЂ4discriminator/discriminator_l2/MatMul/ReadVariableOpЂ5discriminator/discriminator_l3/BiasAdd/ReadVariableOpЂ4discriminator/discriminator_l3/MatMul/ReadVariableOpЂ9discriminator/discriminator_output/BiasAdd/ReadVariableOpЂ8discriminator/discriminator_output/MatMul/ReadVariableOpъ
4discriminator/discriminator_l1/MatMul/ReadVariableOpReadVariableOp=discriminator_discriminator_l1_matmul_readvariableop_resource*
_output_shapes

:
*
dtype026
4discriminator/discriminator_l1/MatMul/ReadVariableOpа
%discriminator/discriminator_l1/MatMulMatMulinputs<discriminator/discriminator_l1/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ
2'
%discriminator/discriminator_l1/MatMulщ
5discriminator/discriminator_l1/BiasAdd/ReadVariableOpReadVariableOp>discriminator_discriminator_l1_biasadd_readvariableop_resource*
_output_shapes
:
*
dtype027
5discriminator/discriminator_l1/BiasAdd/ReadVariableOp§
&discriminator/discriminator_l1/BiasAddBiasAdd/discriminator/discriminator_l1/MatMul:product:0=discriminator/discriminator_l1/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ
2(
&discriminator/discriminator_l1/BiasAddЕ
#discriminator/discriminator_l1/ReluRelu/discriminator/discriminator_l1/BiasAdd:output:0*
T0*'
_output_shapes
:џџџџџџџџџ
2%
#discriminator/discriminator_l1/Reluъ
4discriminator/discriminator_l2/MatMul/ReadVariableOpReadVariableOp=discriminator_discriminator_l2_matmul_readvariableop_resource*
_output_shapes

:
*
dtype026
4discriminator/discriminator_l2/MatMul/ReadVariableOpћ
%discriminator/discriminator_l2/MatMulMatMul1discriminator/discriminator_l1/Relu:activations:0<discriminator/discriminator_l2/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2'
%discriminator/discriminator_l2/MatMulщ
5discriminator/discriminator_l2/BiasAdd/ReadVariableOpReadVariableOp>discriminator_discriminator_l2_biasadd_readvariableop_resource*
_output_shapes
:*
dtype027
5discriminator/discriminator_l2/BiasAdd/ReadVariableOp§
&discriminator/discriminator_l2/BiasAddBiasAdd/discriminator/discriminator_l2/MatMul:product:0=discriminator/discriminator_l2/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2(
&discriminator/discriminator_l2/BiasAddЕ
#discriminator/discriminator_l2/ReluRelu/discriminator/discriminator_l2/BiasAdd:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2%
#discriminator/discriminator_l2/Reluъ
4discriminator/discriminator_l3/MatMul/ReadVariableOpReadVariableOp=discriminator_discriminator_l3_matmul_readvariableop_resource*
_output_shapes

:*
dtype026
4discriminator/discriminator_l3/MatMul/ReadVariableOpћ
%discriminator/discriminator_l3/MatMulMatMul1discriminator/discriminator_l2/Relu:activations:0<discriminator/discriminator_l3/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2'
%discriminator/discriminator_l3/MatMulщ
5discriminator/discriminator_l3/BiasAdd/ReadVariableOpReadVariableOp>discriminator_discriminator_l3_biasadd_readvariableop_resource*
_output_shapes
:*
dtype027
5discriminator/discriminator_l3/BiasAdd/ReadVariableOp§
&discriminator/discriminator_l3/BiasAddBiasAdd/discriminator/discriminator_l3/MatMul:product:0=discriminator/discriminator_l3/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2(
&discriminator/discriminator_l3/BiasAddЕ
#discriminator/discriminator_l3/ReluRelu/discriminator/discriminator_l3/BiasAdd:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2%
#discriminator/discriminator_l3/Reluі
8discriminator/discriminator_output/MatMul/ReadVariableOpReadVariableOpAdiscriminator_discriminator_output_matmul_readvariableop_resource*
_output_shapes

:*
dtype02:
8discriminator/discriminator_output/MatMul/ReadVariableOp
)discriminator/discriminator_output/MatMulMatMul1discriminator/discriminator_l3/Relu:activations:0@discriminator/discriminator_output/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2+
)discriminator/discriminator_output/MatMulѕ
9discriminator/discriminator_output/BiasAdd/ReadVariableOpReadVariableOpBdiscriminator_discriminator_output_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02;
9discriminator/discriminator_output/BiasAdd/ReadVariableOp
*discriminator/discriminator_output/BiasAddBiasAdd3discriminator/discriminator_output/MatMul:product:0Adiscriminator/discriminator_output/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2,
*discriminator/discriminator_output/BiasAddЪ
*discriminator/discriminator_output/SigmoidSigmoid3discriminator/discriminator_output/BiasAdd:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2,
*discriminator/discriminator_output/SigmoidЦ
IdentityIdentity.discriminator/discriminator_output/Sigmoid:y:06^discriminator/discriminator_l1/BiasAdd/ReadVariableOp5^discriminator/discriminator_l1/MatMul/ReadVariableOp6^discriminator/discriminator_l2/BiasAdd/ReadVariableOp5^discriminator/discriminator_l2/MatMul/ReadVariableOp6^discriminator/discriminator_l3/BiasAdd/ReadVariableOp5^discriminator/discriminator_l3/MatMul/ReadVariableOp:^discriminator/discriminator_output/BiasAdd/ReadVariableOp9^discriminator/discriminator_output/MatMul/ReadVariableOp*
T0*'
_output_shapes
:џџџџџџџџџ2

Identity"
identityIdentity:output:0*F
_input_shapes5
3:џџџџџџџџџ::::::::2n
5discriminator/discriminator_l1/BiasAdd/ReadVariableOp5discriminator/discriminator_l1/BiasAdd/ReadVariableOp2l
4discriminator/discriminator_l1/MatMul/ReadVariableOp4discriminator/discriminator_l1/MatMul/ReadVariableOp2n
5discriminator/discriminator_l2/BiasAdd/ReadVariableOp5discriminator/discriminator_l2/BiasAdd/ReadVariableOp2l
4discriminator/discriminator_l2/MatMul/ReadVariableOp4discriminator/discriminator_l2/MatMul/ReadVariableOp2n
5discriminator/discriminator_l3/BiasAdd/ReadVariableOp5discriminator/discriminator_l3/BiasAdd/ReadVariableOp2l
4discriminator/discriminator_l3/MatMul/ReadVariableOp4discriminator/discriminator_l3/MatMul/ReadVariableOp2v
9discriminator/discriminator_output/BiasAdd/ReadVariableOp9discriminator/discriminator_output/BiasAdd/ReadVariableOp2t
8discriminator/discriminator_output/MatMul/ReadVariableOp8discriminator/discriminator_output/MatMul/ReadVariableOp:O K
'
_output_shapes
:џџџџџџџџџ
 
_user_specified_nameinputs

з
'__inference_cgan_31_layer_call_fn_71466
input_1
unknown
	unknown_0
	unknown_1
	unknown_2
	unknown_3
	unknown_4
	unknown_5
	unknown_6
identityЂStatefulPartitionedCallС
StatefulPartitionedCallStatefulPartitionedCallinput_1unknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6*
Tin
2	*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:џџџџџџџџџ**
_read_only_resource_inputs

*-
config_proto

CPU

GPU 2J 8 *K
fFRD
B__inference_cgan_31_layer_call_and_return_conditional_losses_712892
StatefulPartitionedCall
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*'
_output_shapes
:џџџџџџџџџ2

Identity"
identityIdentity:output:0*F
_input_shapes5
3:џџџџџџџџџ::::::::22
StatefulPartitionedCallStatefulPartitionedCall:P L
'
_output_shapes
:џџџџџџџџџ
!
_user_specified_name	input_1

Ы
D__inference_generator_layer_call_and_return_conditional_losses_70888

inputs
generator_l1_70867
generator_l1_70869
generator_l2_70872
generator_l2_70874
generator_l3_70877
generator_l3_70879
generator_output_70882
generator_output_70884
identityЂ$generator_l1/StatefulPartitionedCallЂ$generator_l2/StatefulPartitionedCallЂ$generator_l3/StatefulPartitionedCallЂ(generator_output/StatefulPartitionedCallЅ
$generator_l1/StatefulPartitionedCallStatefulPartitionedCallinputsgenerator_l1_70867generator_l1_70869*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:џџџџџџџџџ
*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8 *P
fKRI
G__inference_generator_l1_layer_call_and_return_conditional_losses_707392&
$generator_l1/StatefulPartitionedCallЬ
$generator_l2/StatefulPartitionedCallStatefulPartitionedCall-generator_l1/StatefulPartitionedCall:output:0generator_l2_70872generator_l2_70874*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:џџџџџџџџџ*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8 *P
fKRI
G__inference_generator_l2_layer_call_and_return_conditional_losses_707662&
$generator_l2/StatefulPartitionedCallЬ
$generator_l3/StatefulPartitionedCallStatefulPartitionedCall-generator_l2/StatefulPartitionedCall:output:0generator_l3_70877generator_l3_70879*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:џџџџџџџџџ*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8 *P
fKRI
G__inference_generator_l3_layer_call_and_return_conditional_losses_707932&
$generator_l3/StatefulPartitionedCallр
(generator_output/StatefulPartitionedCallStatefulPartitionedCall-generator_l3/StatefulPartitionedCall:output:0generator_output_70882generator_output_70884*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:џџџџџџџџџ*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8 *T
fORM
K__inference_generator_output_layer_call_and_return_conditional_losses_708202*
(generator_output/StatefulPartitionedCallЅ
IdentityIdentity1generator_output/StatefulPartitionedCall:output:0%^generator_l1/StatefulPartitionedCall%^generator_l2/StatefulPartitionedCall%^generator_l3/StatefulPartitionedCall)^generator_output/StatefulPartitionedCall*
T0*'
_output_shapes
:џџџџџџџџџ2

Identity"
identityIdentity:output:0*F
_input_shapes5
3:џџџџџџџџџ::::::::2L
$generator_l1/StatefulPartitionedCall$generator_l1/StatefulPartitionedCall2L
$generator_l2/StatefulPartitionedCall$generator_l2/StatefulPartitionedCall2L
$generator_l3/StatefulPartitionedCall$generator_l3/StatefulPartitionedCall2T
(generator_output/StatefulPartitionedCall(generator_output/StatefulPartitionedCall:O K
'
_output_shapes
:џџџџџџџџџ
 
_user_specified_nameinputs
ы

0__inference_discriminator_l1_layer_call_fn_71905

inputs
unknown
	unknown_0
identityЂStatefulPartitionedCallћ
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:џџџџџџџџџ
*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8 *T
fORM
K__inference_discriminator_l1_layer_call_and_return_conditional_losses_709672
StatefulPartitionedCall
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*'
_output_shapes
:џџџџџџџџџ
2

Identity"
identityIdentity:output:0*.
_input_shapes
:џџџџџџџџџ::22
StatefulPartitionedCallStatefulPartitionedCall:O K
'
_output_shapes
:џџџџџџџџџ
 
_user_specified_nameinputs
ы

0__inference_discriminator_l2_layer_call_fn_71925

inputs
unknown
	unknown_0
identityЂStatefulPartitionedCallћ
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:џџџџџџџџџ*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8 *T
fORM
K__inference_discriminator_l2_layer_call_and_return_conditional_losses_709942
StatefulPartitionedCall
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*'
_output_shapes
:џџџџџџџџџ2

Identity"
identityIdentity:output:0*.
_input_shapes
:џџџџџџџџџ
::22
StatefulPartitionedCallStatefulPartitionedCall:O K
'
_output_shapes
:џџџџџџџџџ

 
_user_specified_nameinputs
Ь
щ
-__inference_discriminator_layer_call_fn_71180
discriminator_input
unknown
	unknown_0
	unknown_1
	unknown_2
	unknown_3
	unknown_4
	unknown_5
	unknown_6
identityЂStatefulPartitionedCallг
StatefulPartitionedCallStatefulPartitionedCalldiscriminator_inputunknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6*
Tin
2	*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:џџџџџџџџџ**
_read_only_resource_inputs

*-
config_proto

CPU

GPU 2J 8 *Q
fLRJ
H__inference_discriminator_layer_call_and_return_conditional_losses_711612
StatefulPartitionedCall
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*'
_output_shapes
:џџџџџџџџџ2

Identity"
identityIdentity:output:0*F
_input_shapes5
3:џџџџџџџџџ::::::::22
StatefulPartitionedCallStatefulPartitionedCall:\ X
'
_output_shapes
:џџџџџџџџџ
-
_user_specified_namediscriminator_input
Џ:
Щ
B__inference_cgan_31_layer_call_and_return_conditional_losses_71551

inputsA
=discriminator_discriminator_l1_matmul_readvariableop_resourceB
>discriminator_discriminator_l1_biasadd_readvariableop_resourceA
=discriminator_discriminator_l2_matmul_readvariableop_resourceB
>discriminator_discriminator_l2_biasadd_readvariableop_resourceA
=discriminator_discriminator_l3_matmul_readvariableop_resourceB
>discriminator_discriminator_l3_biasadd_readvariableop_resourceE
Adiscriminator_discriminator_output_matmul_readvariableop_resourceF
Bdiscriminator_discriminator_output_biasadd_readvariableop_resource
identityЂ5discriminator/discriminator_l1/BiasAdd/ReadVariableOpЂ4discriminator/discriminator_l1/MatMul/ReadVariableOpЂ5discriminator/discriminator_l2/BiasAdd/ReadVariableOpЂ4discriminator/discriminator_l2/MatMul/ReadVariableOpЂ5discriminator/discriminator_l3/BiasAdd/ReadVariableOpЂ4discriminator/discriminator_l3/MatMul/ReadVariableOpЂ9discriminator/discriminator_output/BiasAdd/ReadVariableOpЂ8discriminator/discriminator_output/MatMul/ReadVariableOpъ
4discriminator/discriminator_l1/MatMul/ReadVariableOpReadVariableOp=discriminator_discriminator_l1_matmul_readvariableop_resource*
_output_shapes

:
*
dtype026
4discriminator/discriminator_l1/MatMul/ReadVariableOpа
%discriminator/discriminator_l1/MatMulMatMulinputs<discriminator/discriminator_l1/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ
2'
%discriminator/discriminator_l1/MatMulщ
5discriminator/discriminator_l1/BiasAdd/ReadVariableOpReadVariableOp>discriminator_discriminator_l1_biasadd_readvariableop_resource*
_output_shapes
:
*
dtype027
5discriminator/discriminator_l1/BiasAdd/ReadVariableOp§
&discriminator/discriminator_l1/BiasAddBiasAdd/discriminator/discriminator_l1/MatMul:product:0=discriminator/discriminator_l1/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ
2(
&discriminator/discriminator_l1/BiasAddЕ
#discriminator/discriminator_l1/ReluRelu/discriminator/discriminator_l1/BiasAdd:output:0*
T0*'
_output_shapes
:џџџџџџџџџ
2%
#discriminator/discriminator_l1/Reluъ
4discriminator/discriminator_l2/MatMul/ReadVariableOpReadVariableOp=discriminator_discriminator_l2_matmul_readvariableop_resource*
_output_shapes

:
*
dtype026
4discriminator/discriminator_l2/MatMul/ReadVariableOpћ
%discriminator/discriminator_l2/MatMulMatMul1discriminator/discriminator_l1/Relu:activations:0<discriminator/discriminator_l2/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2'
%discriminator/discriminator_l2/MatMulщ
5discriminator/discriminator_l2/BiasAdd/ReadVariableOpReadVariableOp>discriminator_discriminator_l2_biasadd_readvariableop_resource*
_output_shapes
:*
dtype027
5discriminator/discriminator_l2/BiasAdd/ReadVariableOp§
&discriminator/discriminator_l2/BiasAddBiasAdd/discriminator/discriminator_l2/MatMul:product:0=discriminator/discriminator_l2/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2(
&discriminator/discriminator_l2/BiasAddЕ
#discriminator/discriminator_l2/ReluRelu/discriminator/discriminator_l2/BiasAdd:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2%
#discriminator/discriminator_l2/Reluъ
4discriminator/discriminator_l3/MatMul/ReadVariableOpReadVariableOp=discriminator_discriminator_l3_matmul_readvariableop_resource*
_output_shapes

:*
dtype026
4discriminator/discriminator_l3/MatMul/ReadVariableOpћ
%discriminator/discriminator_l3/MatMulMatMul1discriminator/discriminator_l2/Relu:activations:0<discriminator/discriminator_l3/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2'
%discriminator/discriminator_l3/MatMulщ
5discriminator/discriminator_l3/BiasAdd/ReadVariableOpReadVariableOp>discriminator_discriminator_l3_biasadd_readvariableop_resource*
_output_shapes
:*
dtype027
5discriminator/discriminator_l3/BiasAdd/ReadVariableOp§
&discriminator/discriminator_l3/BiasAddBiasAdd/discriminator/discriminator_l3/MatMul:product:0=discriminator/discriminator_l3/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2(
&discriminator/discriminator_l3/BiasAddЕ
#discriminator/discriminator_l3/ReluRelu/discriminator/discriminator_l3/BiasAdd:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2%
#discriminator/discriminator_l3/Reluі
8discriminator/discriminator_output/MatMul/ReadVariableOpReadVariableOpAdiscriminator_discriminator_output_matmul_readvariableop_resource*
_output_shapes

:*
dtype02:
8discriminator/discriminator_output/MatMul/ReadVariableOp
)discriminator/discriminator_output/MatMulMatMul1discriminator/discriminator_l3/Relu:activations:0@discriminator/discriminator_output/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2+
)discriminator/discriminator_output/MatMulѕ
9discriminator/discriminator_output/BiasAdd/ReadVariableOpReadVariableOpBdiscriminator_discriminator_output_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02;
9discriminator/discriminator_output/BiasAdd/ReadVariableOp
*discriminator/discriminator_output/BiasAddBiasAdd3discriminator/discriminator_output/MatMul:product:0Adiscriminator/discriminator_output/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2,
*discriminator/discriminator_output/BiasAddЪ
*discriminator/discriminator_output/SigmoidSigmoid3discriminator/discriminator_output/BiasAdd:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2,
*discriminator/discriminator_output/SigmoidЦ
IdentityIdentity.discriminator/discriminator_output/Sigmoid:y:06^discriminator/discriminator_l1/BiasAdd/ReadVariableOp5^discriminator/discriminator_l1/MatMul/ReadVariableOp6^discriminator/discriminator_l2/BiasAdd/ReadVariableOp5^discriminator/discriminator_l2/MatMul/ReadVariableOp6^discriminator/discriminator_l3/BiasAdd/ReadVariableOp5^discriminator/discriminator_l3/MatMul/ReadVariableOp:^discriminator/discriminator_output/BiasAdd/ReadVariableOp9^discriminator/discriminator_output/MatMul/ReadVariableOp*
T0*'
_output_shapes
:џџџџџџџџџ2

Identity"
identityIdentity:output:0*F
_input_shapes5
3:џџџџџџџџџ::::::::2n
5discriminator/discriminator_l1/BiasAdd/ReadVariableOp5discriminator/discriminator_l1/BiasAdd/ReadVariableOp2l
4discriminator/discriminator_l1/MatMul/ReadVariableOp4discriminator/discriminator_l1/MatMul/ReadVariableOp2n
5discriminator/discriminator_l2/BiasAdd/ReadVariableOp5discriminator/discriminator_l2/BiasAdd/ReadVariableOp2l
4discriminator/discriminator_l2/MatMul/ReadVariableOp4discriminator/discriminator_l2/MatMul/ReadVariableOp2n
5discriminator/discriminator_l3/BiasAdd/ReadVariableOp5discriminator/discriminator_l3/BiasAdd/ReadVariableOp2l
4discriminator/discriminator_l3/MatMul/ReadVariableOp4discriminator/discriminator_l3/MatMul/ReadVariableOp2v
9discriminator/discriminator_output/BiasAdd/ReadVariableOp9discriminator/discriminator_output/BiasAdd/ReadVariableOp2t
8discriminator/discriminator_output/MatMul/ReadVariableOp8discriminator/discriminator_output/MatMul/ReadVariableOp:O K
'
_output_shapes
:џџџџџџџџџ
 
_user_specified_nameinputs

ж
'__inference_cgan_31_layer_call_fn_71572

inputs
unknown
	unknown_0
	unknown_1
	unknown_2
	unknown_3
	unknown_4
	unknown_5
	unknown_6
identityЂStatefulPartitionedCallР
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6*
Tin
2	*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:џџџџџџџџџ**
_read_only_resource_inputs

*-
config_proto

CPU

GPU 2J 8 *K
fFRD
B__inference_cgan_31_layer_call_and_return_conditional_losses_712892
StatefulPartitionedCall
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*'
_output_shapes
:џџџџџџџџџ2

Identity"
identityIdentity:output:0*F
_input_shapes5
3:џџџџџџџџџ::::::::22
StatefulPartitionedCallStatefulPartitionedCall:O K
'
_output_shapes
:џџџџџџџџџ
 
_user_specified_nameinputs
 r
З
__inference__traced_save_72153
file_prefix(
$savev2_adam_iter_read_readvariableop	*
&savev2_adam_beta_1_read_readvariableop*
&savev2_adam_beta_2_read_readvariableop)
%savev2_adam_decay_read_readvariableop1
-savev2_adam_learning_rate_read_readvariableop2
.savev2_generator_l1_kernel_read_readvariableop0
,savev2_generator_l1_bias_read_readvariableop2
.savev2_generator_l2_kernel_read_readvariableop0
,savev2_generator_l2_bias_read_readvariableop2
.savev2_generator_l3_kernel_read_readvariableop0
,savev2_generator_l3_bias_read_readvariableop6
2savev2_generator_output_kernel_read_readvariableop4
0savev2_generator_output_bias_read_readvariableop6
2savev2_discriminator_l1_kernel_read_readvariableop4
0savev2_discriminator_l1_bias_read_readvariableop6
2savev2_discriminator_l2_kernel_read_readvariableop4
0savev2_discriminator_l2_bias_read_readvariableop6
2savev2_discriminator_l3_kernel_read_readvariableop4
0savev2_discriminator_l3_bias_read_readvariableop:
6savev2_discriminator_output_kernel_read_readvariableop8
4savev2_discriminator_output_bias_read_readvariableop$
 savev2_total_read_readvariableop$
 savev2_count_read_readvariableop9
5savev2_adam_generator_l1_kernel_m_read_readvariableop7
3savev2_adam_generator_l1_bias_m_read_readvariableop9
5savev2_adam_generator_l2_kernel_m_read_readvariableop7
3savev2_adam_generator_l2_bias_m_read_readvariableop9
5savev2_adam_generator_l3_kernel_m_read_readvariableop7
3savev2_adam_generator_l3_bias_m_read_readvariableop=
9savev2_adam_generator_output_kernel_m_read_readvariableop;
7savev2_adam_generator_output_bias_m_read_readvariableop?
;savev2_adam_1_discriminator_l1_kernel_m_read_readvariableop=
9savev2_adam_1_discriminator_l1_bias_m_read_readvariableop?
;savev2_adam_1_discriminator_l2_kernel_m_read_readvariableop=
9savev2_adam_1_discriminator_l2_bias_m_read_readvariableop?
;savev2_adam_1_discriminator_l3_kernel_m_read_readvariableop=
9savev2_adam_1_discriminator_l3_bias_m_read_readvariableopC
?savev2_adam_1_discriminator_output_kernel_m_read_readvariableopA
=savev2_adam_1_discriminator_output_bias_m_read_readvariableop9
5savev2_adam_generator_l1_kernel_v_read_readvariableop7
3savev2_adam_generator_l1_bias_v_read_readvariableop9
5savev2_adam_generator_l2_kernel_v_read_readvariableop7
3savev2_adam_generator_l2_bias_v_read_readvariableop9
5savev2_adam_generator_l3_kernel_v_read_readvariableop7
3savev2_adam_generator_l3_bias_v_read_readvariableop=
9savev2_adam_generator_output_kernel_v_read_readvariableop;
7savev2_adam_generator_output_bias_v_read_readvariableop?
;savev2_adam_1_discriminator_l1_kernel_v_read_readvariableop=
9savev2_adam_1_discriminator_l1_bias_v_read_readvariableop?
;savev2_adam_1_discriminator_l2_kernel_v_read_readvariableop=
9savev2_adam_1_discriminator_l2_bias_v_read_readvariableop?
;savev2_adam_1_discriminator_l3_kernel_v_read_readvariableop=
9savev2_adam_1_discriminator_l3_bias_v_read_readvariableopC
?savev2_adam_1_discriminator_output_kernel_v_read_readvariableopA
=savev2_adam_1_discriminator_output_bias_v_read_readvariableop
savev2_const

identity_1ЂMergeV2Checkpoints
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
Const_1
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
ShardedFilename/shardІ
ShardedFilenameShardedFilenameStringJoin:output:0ShardedFilename/shard:output:0num_shards:output:0"/device:CPU:0*
_output_shapes
: 2
ShardedFilenameф
SaveV2/tensor_namesConst"/device:CPU:0*
_output_shapes
:8*
dtype0*і
valueьBщ8B)optimizer/iter/.ATTRIBUTES/VARIABLE_VALUEB+optimizer/beta_1/.ATTRIBUTES/VARIABLE_VALUEB+optimizer/beta_2/.ATTRIBUTES/VARIABLE_VALUEB*optimizer/decay/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/learning_rate/.ATTRIBUTES/VARIABLE_VALUEB0trainable_variables/0/.ATTRIBUTES/VARIABLE_VALUEB0trainable_variables/1/.ATTRIBUTES/VARIABLE_VALUEB0trainable_variables/2/.ATTRIBUTES/VARIABLE_VALUEB0trainable_variables/3/.ATTRIBUTES/VARIABLE_VALUEB0trainable_variables/4/.ATTRIBUTES/VARIABLE_VALUEB0trainable_variables/5/.ATTRIBUTES/VARIABLE_VALUEB0trainable_variables/6/.ATTRIBUTES/VARIABLE_VALUEB0trainable_variables/7/.ATTRIBUTES/VARIABLE_VALUEB0trainable_variables/8/.ATTRIBUTES/VARIABLE_VALUEB0trainable_variables/9/.ATTRIBUTES/VARIABLE_VALUEB1trainable_variables/10/.ATTRIBUTES/VARIABLE_VALUEB1trainable_variables/11/.ATTRIBUTES/VARIABLE_VALUEB1trainable_variables/12/.ATTRIBUTES/VARIABLE_VALUEB1trainable_variables/13/.ATTRIBUTES/VARIABLE_VALUEB1trainable_variables/14/.ATTRIBUTES/VARIABLE_VALUEB1trainable_variables/15/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/0/total/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/0/count/.ATTRIBUTES/VARIABLE_VALUEBLtrainable_variables/0/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBLtrainable_variables/1/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBLtrainable_variables/2/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBLtrainable_variables/3/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBLtrainable_variables/4/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBLtrainable_variables/5/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBLtrainable_variables/6/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBLtrainable_variables/7/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBLtrainable_variables/8/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBLtrainable_variables/9/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBMtrainable_variables/10/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBMtrainable_variables/11/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBMtrainable_variables/12/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBMtrainable_variables/13/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBMtrainable_variables/14/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBMtrainable_variables/15/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBLtrainable_variables/0/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBLtrainable_variables/1/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBLtrainable_variables/2/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBLtrainable_variables/3/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBLtrainable_variables/4/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBLtrainable_variables/5/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBLtrainable_variables/6/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBLtrainable_variables/7/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBLtrainable_variables/8/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBLtrainable_variables/9/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBMtrainable_variables/10/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBMtrainable_variables/11/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBMtrainable_variables/12/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBMtrainable_variables/13/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBMtrainable_variables/14/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBMtrainable_variables/15/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEB_CHECKPOINTABLE_OBJECT_GRAPH2
SaveV2/tensor_namesљ
SaveV2/shape_and_slicesConst"/device:CPU:0*
_output_shapes
:8*
dtype0*
valuezBx8B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B 2
SaveV2/shape_and_slicesс
SaveV2SaveV2ShardedFilename:filename:0SaveV2/tensor_names:output:0 SaveV2/shape_and_slices:output:0$savev2_adam_iter_read_readvariableop&savev2_adam_beta_1_read_readvariableop&savev2_adam_beta_2_read_readvariableop%savev2_adam_decay_read_readvariableop-savev2_adam_learning_rate_read_readvariableop.savev2_generator_l1_kernel_read_readvariableop,savev2_generator_l1_bias_read_readvariableop.savev2_generator_l2_kernel_read_readvariableop,savev2_generator_l2_bias_read_readvariableop.savev2_generator_l3_kernel_read_readvariableop,savev2_generator_l3_bias_read_readvariableop2savev2_generator_output_kernel_read_readvariableop0savev2_generator_output_bias_read_readvariableop2savev2_discriminator_l1_kernel_read_readvariableop0savev2_discriminator_l1_bias_read_readvariableop2savev2_discriminator_l2_kernel_read_readvariableop0savev2_discriminator_l2_bias_read_readvariableop2savev2_discriminator_l3_kernel_read_readvariableop0savev2_discriminator_l3_bias_read_readvariableop6savev2_discriminator_output_kernel_read_readvariableop4savev2_discriminator_output_bias_read_readvariableop savev2_total_read_readvariableop savev2_count_read_readvariableop5savev2_adam_generator_l1_kernel_m_read_readvariableop3savev2_adam_generator_l1_bias_m_read_readvariableop5savev2_adam_generator_l2_kernel_m_read_readvariableop3savev2_adam_generator_l2_bias_m_read_readvariableop5savev2_adam_generator_l3_kernel_m_read_readvariableop3savev2_adam_generator_l3_bias_m_read_readvariableop9savev2_adam_generator_output_kernel_m_read_readvariableop7savev2_adam_generator_output_bias_m_read_readvariableop;savev2_adam_1_discriminator_l1_kernel_m_read_readvariableop9savev2_adam_1_discriminator_l1_bias_m_read_readvariableop;savev2_adam_1_discriminator_l2_kernel_m_read_readvariableop9savev2_adam_1_discriminator_l2_bias_m_read_readvariableop;savev2_adam_1_discriminator_l3_kernel_m_read_readvariableop9savev2_adam_1_discriminator_l3_bias_m_read_readvariableop?savev2_adam_1_discriminator_output_kernel_m_read_readvariableop=savev2_adam_1_discriminator_output_bias_m_read_readvariableop5savev2_adam_generator_l1_kernel_v_read_readvariableop3savev2_adam_generator_l1_bias_v_read_readvariableop5savev2_adam_generator_l2_kernel_v_read_readvariableop3savev2_adam_generator_l2_bias_v_read_readvariableop5savev2_adam_generator_l3_kernel_v_read_readvariableop3savev2_adam_generator_l3_bias_v_read_readvariableop9savev2_adam_generator_output_kernel_v_read_readvariableop7savev2_adam_generator_output_bias_v_read_readvariableop;savev2_adam_1_discriminator_l1_kernel_v_read_readvariableop9savev2_adam_1_discriminator_l1_bias_v_read_readvariableop;savev2_adam_1_discriminator_l2_kernel_v_read_readvariableop9savev2_adam_1_discriminator_l2_bias_v_read_readvariableop;savev2_adam_1_discriminator_l3_kernel_v_read_readvariableop9savev2_adam_1_discriminator_l3_bias_v_read_readvariableop?savev2_adam_1_discriminator_output_kernel_v_read_readvariableop=savev2_adam_1_discriminator_output_bias_v_read_readvariableopsavev2_const"/device:CPU:0*
_output_shapes
 *F
dtypes<
:28	2
SaveV2К
&MergeV2Checkpoints/checkpoint_prefixesPackShardedFilename:filename:0^SaveV2"/device:CPU:0*
N*
T0*
_output_shapes
:2(
&MergeV2Checkpoints/checkpoint_prefixesЁ
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

identity_1Identity_1:output:0*Ї
_input_shapes
: : : : : : :
:
:
::::::
:
:
:::::: : :
:
:
::::::
:
:
::::::
:
:
::::::
:
:
:::::: 2(
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

:
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

:: 

_output_shapes
::$ 

_output_shapes

:
: 

_output_shapes
:
:$ 

_output_shapes

:
: 

_output_shapes
::$ 

_output_shapes

:: 

_output_shapes
::$ 

_output_shapes

:: 

_output_shapes
::

_output_shapes
: :

_output_shapes
: :$ 

_output_shapes

:
: 

_output_shapes
:
:$ 

_output_shapes

:
: 

_output_shapes
::$ 

_output_shapes

:: 

_output_shapes
::$ 

_output_shapes

:: 

_output_shapes
::$  

_output_shapes

:
: !

_output_shapes
:
:$" 

_output_shapes

:
: #

_output_shapes
::$$ 

_output_shapes

:: %

_output_shapes
::$& 

_output_shapes

:: '

_output_shapes
::$( 

_output_shapes

:
: )

_output_shapes
:
:$* 

_output_shapes

:
: +

_output_shapes
::$, 

_output_shapes

:: -

_output_shapes
::$. 

_output_shapes

:: /

_output_shapes
::$0 

_output_shapes

:
: 1

_output_shapes
:
:$2 

_output_shapes

:
: 3

_output_shapes
::$4 

_output_shapes

:: 5

_output_shapes
::$6 

_output_shapes

:: 7

_output_shapes
::8

_output_shapes
: 
у

,__inference_generator_l2_layer_call_fn_71845

inputs
unknown
	unknown_0
identityЂStatefulPartitionedCallї
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:џџџџџџџџџ*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8 *P
fKRI
G__inference_generator_l2_layer_call_and_return_conditional_losses_707662
StatefulPartitionedCall
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*'
_output_shapes
:џџџџџџџџџ2

Identity"
identityIdentity:output:0*.
_input_shapes
:џџџџџџџџџ
::22
StatefulPartitionedCallStatefulPartitionedCall:O K
'
_output_shapes
:џџџџџџџџџ

 
_user_specified_nameinputs

и
)__inference_generator_layer_call_fn_71699

inputs
unknown
	unknown_0
	unknown_1
	unknown_2
	unknown_3
	unknown_4
	unknown_5
	unknown_6
identityЂStatefulPartitionedCallТ
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6*
Tin
2	*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:џџџџџџџџџ**
_read_only_resource_inputs

*-
config_proto

CPU

GPU 2J 8 *M
fHRF
D__inference_generator_layer_call_and_return_conditional_losses_709332
StatefulPartitionedCall
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*'
_output_shapes
:џџџџџџџџџ2

Identity"
identityIdentity:output:0*F
_input_shapes5
3:џџџџџџџџџ::::::::22
StatefulPartitionedCallStatefulPartitionedCall:O K
'
_output_shapes
:џџџџџџџџџ
 
_user_specified_nameinputs
Ќы
Е
!__inference__traced_restore_72328
file_prefix
assignvariableop_adam_iter"
assignvariableop_1_adam_beta_1"
assignvariableop_2_adam_beta_2!
assignvariableop_3_adam_decay)
%assignvariableop_4_adam_learning_rate*
&assignvariableop_5_generator_l1_kernel(
$assignvariableop_6_generator_l1_bias*
&assignvariableop_7_generator_l2_kernel(
$assignvariableop_8_generator_l2_bias*
&assignvariableop_9_generator_l3_kernel)
%assignvariableop_10_generator_l3_bias/
+assignvariableop_11_generator_output_kernel-
)assignvariableop_12_generator_output_bias/
+assignvariableop_13_discriminator_l1_kernel-
)assignvariableop_14_discriminator_l1_bias/
+assignvariableop_15_discriminator_l2_kernel-
)assignvariableop_16_discriminator_l2_bias/
+assignvariableop_17_discriminator_l3_kernel-
)assignvariableop_18_discriminator_l3_bias3
/assignvariableop_19_discriminator_output_kernel1
-assignvariableop_20_discriminator_output_bias
assignvariableop_21_total
assignvariableop_22_count2
.assignvariableop_23_adam_generator_l1_kernel_m0
,assignvariableop_24_adam_generator_l1_bias_m2
.assignvariableop_25_adam_generator_l2_kernel_m0
,assignvariableop_26_adam_generator_l2_bias_m2
.assignvariableop_27_adam_generator_l3_kernel_m0
,assignvariableop_28_adam_generator_l3_bias_m6
2assignvariableop_29_adam_generator_output_kernel_m4
0assignvariableop_30_adam_generator_output_bias_m8
4assignvariableop_31_adam_1_discriminator_l1_kernel_m6
2assignvariableop_32_adam_1_discriminator_l1_bias_m8
4assignvariableop_33_adam_1_discriminator_l2_kernel_m6
2assignvariableop_34_adam_1_discriminator_l2_bias_m8
4assignvariableop_35_adam_1_discriminator_l3_kernel_m6
2assignvariableop_36_adam_1_discriminator_l3_bias_m<
8assignvariableop_37_adam_1_discriminator_output_kernel_m:
6assignvariableop_38_adam_1_discriminator_output_bias_m2
.assignvariableop_39_adam_generator_l1_kernel_v0
,assignvariableop_40_adam_generator_l1_bias_v2
.assignvariableop_41_adam_generator_l2_kernel_v0
,assignvariableop_42_adam_generator_l2_bias_v2
.assignvariableop_43_adam_generator_l3_kernel_v0
,assignvariableop_44_adam_generator_l3_bias_v6
2assignvariableop_45_adam_generator_output_kernel_v4
0assignvariableop_46_adam_generator_output_bias_v8
4assignvariableop_47_adam_1_discriminator_l1_kernel_v6
2assignvariableop_48_adam_1_discriminator_l1_bias_v8
4assignvariableop_49_adam_1_discriminator_l2_kernel_v6
2assignvariableop_50_adam_1_discriminator_l2_bias_v8
4assignvariableop_51_adam_1_discriminator_l3_kernel_v6
2assignvariableop_52_adam_1_discriminator_l3_bias_v<
8assignvariableop_53_adam_1_discriminator_output_kernel_v:
6assignvariableop_54_adam_1_discriminator_output_bias_v
identity_56ЂAssignVariableOpЂAssignVariableOp_1ЂAssignVariableOp_10ЂAssignVariableOp_11ЂAssignVariableOp_12ЂAssignVariableOp_13ЂAssignVariableOp_14ЂAssignVariableOp_15ЂAssignVariableOp_16ЂAssignVariableOp_17ЂAssignVariableOp_18ЂAssignVariableOp_19ЂAssignVariableOp_2ЂAssignVariableOp_20ЂAssignVariableOp_21ЂAssignVariableOp_22ЂAssignVariableOp_23ЂAssignVariableOp_24ЂAssignVariableOp_25ЂAssignVariableOp_26ЂAssignVariableOp_27ЂAssignVariableOp_28ЂAssignVariableOp_29ЂAssignVariableOp_3ЂAssignVariableOp_30ЂAssignVariableOp_31ЂAssignVariableOp_32ЂAssignVariableOp_33ЂAssignVariableOp_34ЂAssignVariableOp_35ЂAssignVariableOp_36ЂAssignVariableOp_37ЂAssignVariableOp_38ЂAssignVariableOp_39ЂAssignVariableOp_4ЂAssignVariableOp_40ЂAssignVariableOp_41ЂAssignVariableOp_42ЂAssignVariableOp_43ЂAssignVariableOp_44ЂAssignVariableOp_45ЂAssignVariableOp_46ЂAssignVariableOp_47ЂAssignVariableOp_48ЂAssignVariableOp_49ЂAssignVariableOp_5ЂAssignVariableOp_50ЂAssignVariableOp_51ЂAssignVariableOp_52ЂAssignVariableOp_53ЂAssignVariableOp_54ЂAssignVariableOp_6ЂAssignVariableOp_7ЂAssignVariableOp_8ЂAssignVariableOp_9ъ
RestoreV2/tensor_namesConst"/device:CPU:0*
_output_shapes
:8*
dtype0*і
valueьBщ8B)optimizer/iter/.ATTRIBUTES/VARIABLE_VALUEB+optimizer/beta_1/.ATTRIBUTES/VARIABLE_VALUEB+optimizer/beta_2/.ATTRIBUTES/VARIABLE_VALUEB*optimizer/decay/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/learning_rate/.ATTRIBUTES/VARIABLE_VALUEB0trainable_variables/0/.ATTRIBUTES/VARIABLE_VALUEB0trainable_variables/1/.ATTRIBUTES/VARIABLE_VALUEB0trainable_variables/2/.ATTRIBUTES/VARIABLE_VALUEB0trainable_variables/3/.ATTRIBUTES/VARIABLE_VALUEB0trainable_variables/4/.ATTRIBUTES/VARIABLE_VALUEB0trainable_variables/5/.ATTRIBUTES/VARIABLE_VALUEB0trainable_variables/6/.ATTRIBUTES/VARIABLE_VALUEB0trainable_variables/7/.ATTRIBUTES/VARIABLE_VALUEB0trainable_variables/8/.ATTRIBUTES/VARIABLE_VALUEB0trainable_variables/9/.ATTRIBUTES/VARIABLE_VALUEB1trainable_variables/10/.ATTRIBUTES/VARIABLE_VALUEB1trainable_variables/11/.ATTRIBUTES/VARIABLE_VALUEB1trainable_variables/12/.ATTRIBUTES/VARIABLE_VALUEB1trainable_variables/13/.ATTRIBUTES/VARIABLE_VALUEB1trainable_variables/14/.ATTRIBUTES/VARIABLE_VALUEB1trainable_variables/15/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/0/total/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/0/count/.ATTRIBUTES/VARIABLE_VALUEBLtrainable_variables/0/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBLtrainable_variables/1/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBLtrainable_variables/2/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBLtrainable_variables/3/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBLtrainable_variables/4/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBLtrainable_variables/5/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBLtrainable_variables/6/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBLtrainable_variables/7/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBLtrainable_variables/8/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBLtrainable_variables/9/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBMtrainable_variables/10/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBMtrainable_variables/11/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBMtrainable_variables/12/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBMtrainable_variables/13/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBMtrainable_variables/14/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBMtrainable_variables/15/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBLtrainable_variables/0/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBLtrainable_variables/1/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBLtrainable_variables/2/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBLtrainable_variables/3/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBLtrainable_variables/4/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBLtrainable_variables/5/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBLtrainable_variables/6/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBLtrainable_variables/7/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBLtrainable_variables/8/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBLtrainable_variables/9/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBMtrainable_variables/10/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBMtrainable_variables/11/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBMtrainable_variables/12/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBMtrainable_variables/13/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBMtrainable_variables/14/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBMtrainable_variables/15/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEB_CHECKPOINTABLE_OBJECT_GRAPH2
RestoreV2/tensor_namesџ
RestoreV2/shape_and_slicesConst"/device:CPU:0*
_output_shapes
:8*
dtype0*
valuezBx8B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B 2
RestoreV2/shape_and_slicesЦ
	RestoreV2	RestoreV2file_prefixRestoreV2/tensor_names:output:0#RestoreV2/shape_and_slices:output:0"/device:CPU:0*і
_output_shapesу
р::::::::::::::::::::::::::::::::::::::::::::::::::::::::*F
dtypes<
:28	2
	RestoreV2g
IdentityIdentityRestoreV2:tensors:0"/device:CPU:0*
T0	*
_output_shapes
:2

Identity
AssignVariableOpAssignVariableOpassignvariableop_adam_iterIdentity:output:0"/device:CPU:0*
_output_shapes
 *
dtype0	2
AssignVariableOpk

Identity_1IdentityRestoreV2:tensors:1"/device:CPU:0*
T0*
_output_shapes
:2

Identity_1Ѓ
AssignVariableOp_1AssignVariableOpassignvariableop_1_adam_beta_1Identity_1:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_1k

Identity_2IdentityRestoreV2:tensors:2"/device:CPU:0*
T0*
_output_shapes
:2

Identity_2Ѓ
AssignVariableOp_2AssignVariableOpassignvariableop_2_adam_beta_2Identity_2:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_2k

Identity_3IdentityRestoreV2:tensors:3"/device:CPU:0*
T0*
_output_shapes
:2

Identity_3Ђ
AssignVariableOp_3AssignVariableOpassignvariableop_3_adam_decayIdentity_3:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_3k

Identity_4IdentityRestoreV2:tensors:4"/device:CPU:0*
T0*
_output_shapes
:2

Identity_4Њ
AssignVariableOp_4AssignVariableOp%assignvariableop_4_adam_learning_rateIdentity_4:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_4k

Identity_5IdentityRestoreV2:tensors:5"/device:CPU:0*
T0*
_output_shapes
:2

Identity_5Ћ
AssignVariableOp_5AssignVariableOp&assignvariableop_5_generator_l1_kernelIdentity_5:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_5k

Identity_6IdentityRestoreV2:tensors:6"/device:CPU:0*
T0*
_output_shapes
:2

Identity_6Љ
AssignVariableOp_6AssignVariableOp$assignvariableop_6_generator_l1_biasIdentity_6:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_6k

Identity_7IdentityRestoreV2:tensors:7"/device:CPU:0*
T0*
_output_shapes
:2

Identity_7Ћ
AssignVariableOp_7AssignVariableOp&assignvariableop_7_generator_l2_kernelIdentity_7:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_7k

Identity_8IdentityRestoreV2:tensors:8"/device:CPU:0*
T0*
_output_shapes
:2

Identity_8Љ
AssignVariableOp_8AssignVariableOp$assignvariableop_8_generator_l2_biasIdentity_8:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_8k

Identity_9IdentityRestoreV2:tensors:9"/device:CPU:0*
T0*
_output_shapes
:2

Identity_9Ћ
AssignVariableOp_9AssignVariableOp&assignvariableop_9_generator_l3_kernelIdentity_9:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_9n
Identity_10IdentityRestoreV2:tensors:10"/device:CPU:0*
T0*
_output_shapes
:2
Identity_10­
AssignVariableOp_10AssignVariableOp%assignvariableop_10_generator_l3_biasIdentity_10:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_10n
Identity_11IdentityRestoreV2:tensors:11"/device:CPU:0*
T0*
_output_shapes
:2
Identity_11Г
AssignVariableOp_11AssignVariableOp+assignvariableop_11_generator_output_kernelIdentity_11:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_11n
Identity_12IdentityRestoreV2:tensors:12"/device:CPU:0*
T0*
_output_shapes
:2
Identity_12Б
AssignVariableOp_12AssignVariableOp)assignvariableop_12_generator_output_biasIdentity_12:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_12n
Identity_13IdentityRestoreV2:tensors:13"/device:CPU:0*
T0*
_output_shapes
:2
Identity_13Г
AssignVariableOp_13AssignVariableOp+assignvariableop_13_discriminator_l1_kernelIdentity_13:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_13n
Identity_14IdentityRestoreV2:tensors:14"/device:CPU:0*
T0*
_output_shapes
:2
Identity_14Б
AssignVariableOp_14AssignVariableOp)assignvariableop_14_discriminator_l1_biasIdentity_14:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_14n
Identity_15IdentityRestoreV2:tensors:15"/device:CPU:0*
T0*
_output_shapes
:2
Identity_15Г
AssignVariableOp_15AssignVariableOp+assignvariableop_15_discriminator_l2_kernelIdentity_15:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_15n
Identity_16IdentityRestoreV2:tensors:16"/device:CPU:0*
T0*
_output_shapes
:2
Identity_16Б
AssignVariableOp_16AssignVariableOp)assignvariableop_16_discriminator_l2_biasIdentity_16:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_16n
Identity_17IdentityRestoreV2:tensors:17"/device:CPU:0*
T0*
_output_shapes
:2
Identity_17Г
AssignVariableOp_17AssignVariableOp+assignvariableop_17_discriminator_l3_kernelIdentity_17:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_17n
Identity_18IdentityRestoreV2:tensors:18"/device:CPU:0*
T0*
_output_shapes
:2
Identity_18Б
AssignVariableOp_18AssignVariableOp)assignvariableop_18_discriminator_l3_biasIdentity_18:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_18n
Identity_19IdentityRestoreV2:tensors:19"/device:CPU:0*
T0*
_output_shapes
:2
Identity_19З
AssignVariableOp_19AssignVariableOp/assignvariableop_19_discriminator_output_kernelIdentity_19:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_19n
Identity_20IdentityRestoreV2:tensors:20"/device:CPU:0*
T0*
_output_shapes
:2
Identity_20Е
AssignVariableOp_20AssignVariableOp-assignvariableop_20_discriminator_output_biasIdentity_20:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_20n
Identity_21IdentityRestoreV2:tensors:21"/device:CPU:0*
T0*
_output_shapes
:2
Identity_21Ё
AssignVariableOp_21AssignVariableOpassignvariableop_21_totalIdentity_21:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_21n
Identity_22IdentityRestoreV2:tensors:22"/device:CPU:0*
T0*
_output_shapes
:2
Identity_22Ё
AssignVariableOp_22AssignVariableOpassignvariableop_22_countIdentity_22:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_22n
Identity_23IdentityRestoreV2:tensors:23"/device:CPU:0*
T0*
_output_shapes
:2
Identity_23Ж
AssignVariableOp_23AssignVariableOp.assignvariableop_23_adam_generator_l1_kernel_mIdentity_23:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_23n
Identity_24IdentityRestoreV2:tensors:24"/device:CPU:0*
T0*
_output_shapes
:2
Identity_24Д
AssignVariableOp_24AssignVariableOp,assignvariableop_24_adam_generator_l1_bias_mIdentity_24:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_24n
Identity_25IdentityRestoreV2:tensors:25"/device:CPU:0*
T0*
_output_shapes
:2
Identity_25Ж
AssignVariableOp_25AssignVariableOp.assignvariableop_25_adam_generator_l2_kernel_mIdentity_25:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_25n
Identity_26IdentityRestoreV2:tensors:26"/device:CPU:0*
T0*
_output_shapes
:2
Identity_26Д
AssignVariableOp_26AssignVariableOp,assignvariableop_26_adam_generator_l2_bias_mIdentity_26:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_26n
Identity_27IdentityRestoreV2:tensors:27"/device:CPU:0*
T0*
_output_shapes
:2
Identity_27Ж
AssignVariableOp_27AssignVariableOp.assignvariableop_27_adam_generator_l3_kernel_mIdentity_27:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_27n
Identity_28IdentityRestoreV2:tensors:28"/device:CPU:0*
T0*
_output_shapes
:2
Identity_28Д
AssignVariableOp_28AssignVariableOp,assignvariableop_28_adam_generator_l3_bias_mIdentity_28:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_28n
Identity_29IdentityRestoreV2:tensors:29"/device:CPU:0*
T0*
_output_shapes
:2
Identity_29К
AssignVariableOp_29AssignVariableOp2assignvariableop_29_adam_generator_output_kernel_mIdentity_29:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_29n
Identity_30IdentityRestoreV2:tensors:30"/device:CPU:0*
T0*
_output_shapes
:2
Identity_30И
AssignVariableOp_30AssignVariableOp0assignvariableop_30_adam_generator_output_bias_mIdentity_30:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_30n
Identity_31IdentityRestoreV2:tensors:31"/device:CPU:0*
T0*
_output_shapes
:2
Identity_31М
AssignVariableOp_31AssignVariableOp4assignvariableop_31_adam_1_discriminator_l1_kernel_mIdentity_31:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_31n
Identity_32IdentityRestoreV2:tensors:32"/device:CPU:0*
T0*
_output_shapes
:2
Identity_32К
AssignVariableOp_32AssignVariableOp2assignvariableop_32_adam_1_discriminator_l1_bias_mIdentity_32:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_32n
Identity_33IdentityRestoreV2:tensors:33"/device:CPU:0*
T0*
_output_shapes
:2
Identity_33М
AssignVariableOp_33AssignVariableOp4assignvariableop_33_adam_1_discriminator_l2_kernel_mIdentity_33:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_33n
Identity_34IdentityRestoreV2:tensors:34"/device:CPU:0*
T0*
_output_shapes
:2
Identity_34К
AssignVariableOp_34AssignVariableOp2assignvariableop_34_adam_1_discriminator_l2_bias_mIdentity_34:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_34n
Identity_35IdentityRestoreV2:tensors:35"/device:CPU:0*
T0*
_output_shapes
:2
Identity_35М
AssignVariableOp_35AssignVariableOp4assignvariableop_35_adam_1_discriminator_l3_kernel_mIdentity_35:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_35n
Identity_36IdentityRestoreV2:tensors:36"/device:CPU:0*
T0*
_output_shapes
:2
Identity_36К
AssignVariableOp_36AssignVariableOp2assignvariableop_36_adam_1_discriminator_l3_bias_mIdentity_36:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_36n
Identity_37IdentityRestoreV2:tensors:37"/device:CPU:0*
T0*
_output_shapes
:2
Identity_37Р
AssignVariableOp_37AssignVariableOp8assignvariableop_37_adam_1_discriminator_output_kernel_mIdentity_37:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_37n
Identity_38IdentityRestoreV2:tensors:38"/device:CPU:0*
T0*
_output_shapes
:2
Identity_38О
AssignVariableOp_38AssignVariableOp6assignvariableop_38_adam_1_discriminator_output_bias_mIdentity_38:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_38n
Identity_39IdentityRestoreV2:tensors:39"/device:CPU:0*
T0*
_output_shapes
:2
Identity_39Ж
AssignVariableOp_39AssignVariableOp.assignvariableop_39_adam_generator_l1_kernel_vIdentity_39:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_39n
Identity_40IdentityRestoreV2:tensors:40"/device:CPU:0*
T0*
_output_shapes
:2
Identity_40Д
AssignVariableOp_40AssignVariableOp,assignvariableop_40_adam_generator_l1_bias_vIdentity_40:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_40n
Identity_41IdentityRestoreV2:tensors:41"/device:CPU:0*
T0*
_output_shapes
:2
Identity_41Ж
AssignVariableOp_41AssignVariableOp.assignvariableop_41_adam_generator_l2_kernel_vIdentity_41:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_41n
Identity_42IdentityRestoreV2:tensors:42"/device:CPU:0*
T0*
_output_shapes
:2
Identity_42Д
AssignVariableOp_42AssignVariableOp,assignvariableop_42_adam_generator_l2_bias_vIdentity_42:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_42n
Identity_43IdentityRestoreV2:tensors:43"/device:CPU:0*
T0*
_output_shapes
:2
Identity_43Ж
AssignVariableOp_43AssignVariableOp.assignvariableop_43_adam_generator_l3_kernel_vIdentity_43:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_43n
Identity_44IdentityRestoreV2:tensors:44"/device:CPU:0*
T0*
_output_shapes
:2
Identity_44Д
AssignVariableOp_44AssignVariableOp,assignvariableop_44_adam_generator_l3_bias_vIdentity_44:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_44n
Identity_45IdentityRestoreV2:tensors:45"/device:CPU:0*
T0*
_output_shapes
:2
Identity_45К
AssignVariableOp_45AssignVariableOp2assignvariableop_45_adam_generator_output_kernel_vIdentity_45:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_45n
Identity_46IdentityRestoreV2:tensors:46"/device:CPU:0*
T0*
_output_shapes
:2
Identity_46И
AssignVariableOp_46AssignVariableOp0assignvariableop_46_adam_generator_output_bias_vIdentity_46:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_46n
Identity_47IdentityRestoreV2:tensors:47"/device:CPU:0*
T0*
_output_shapes
:2
Identity_47М
AssignVariableOp_47AssignVariableOp4assignvariableop_47_adam_1_discriminator_l1_kernel_vIdentity_47:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_47n
Identity_48IdentityRestoreV2:tensors:48"/device:CPU:0*
T0*
_output_shapes
:2
Identity_48К
AssignVariableOp_48AssignVariableOp2assignvariableop_48_adam_1_discriminator_l1_bias_vIdentity_48:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_48n
Identity_49IdentityRestoreV2:tensors:49"/device:CPU:0*
T0*
_output_shapes
:2
Identity_49М
AssignVariableOp_49AssignVariableOp4assignvariableop_49_adam_1_discriminator_l2_kernel_vIdentity_49:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_49n
Identity_50IdentityRestoreV2:tensors:50"/device:CPU:0*
T0*
_output_shapes
:2
Identity_50К
AssignVariableOp_50AssignVariableOp2assignvariableop_50_adam_1_discriminator_l2_bias_vIdentity_50:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_50n
Identity_51IdentityRestoreV2:tensors:51"/device:CPU:0*
T0*
_output_shapes
:2
Identity_51М
AssignVariableOp_51AssignVariableOp4assignvariableop_51_adam_1_discriminator_l3_kernel_vIdentity_51:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_51n
Identity_52IdentityRestoreV2:tensors:52"/device:CPU:0*
T0*
_output_shapes
:2
Identity_52К
AssignVariableOp_52AssignVariableOp2assignvariableop_52_adam_1_discriminator_l3_bias_vIdentity_52:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_52n
Identity_53IdentityRestoreV2:tensors:53"/device:CPU:0*
T0*
_output_shapes
:2
Identity_53Р
AssignVariableOp_53AssignVariableOp8assignvariableop_53_adam_1_discriminator_output_kernel_vIdentity_53:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_53n
Identity_54IdentityRestoreV2:tensors:54"/device:CPU:0*
T0*
_output_shapes
:2
Identity_54О
AssignVariableOp_54AssignVariableOp6assignvariableop_54_adam_1_discriminator_output_bias_vIdentity_54:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_549
NoOpNoOp"/device:CPU:0*
_output_shapes
 2
NoOp

Identity_55Identityfile_prefix^AssignVariableOp^AssignVariableOp_1^AssignVariableOp_10^AssignVariableOp_11^AssignVariableOp_12^AssignVariableOp_13^AssignVariableOp_14^AssignVariableOp_15^AssignVariableOp_16^AssignVariableOp_17^AssignVariableOp_18^AssignVariableOp_19^AssignVariableOp_2^AssignVariableOp_20^AssignVariableOp_21^AssignVariableOp_22^AssignVariableOp_23^AssignVariableOp_24^AssignVariableOp_25^AssignVariableOp_26^AssignVariableOp_27^AssignVariableOp_28^AssignVariableOp_29^AssignVariableOp_3^AssignVariableOp_30^AssignVariableOp_31^AssignVariableOp_32^AssignVariableOp_33^AssignVariableOp_34^AssignVariableOp_35^AssignVariableOp_36^AssignVariableOp_37^AssignVariableOp_38^AssignVariableOp_39^AssignVariableOp_4^AssignVariableOp_40^AssignVariableOp_41^AssignVariableOp_42^AssignVariableOp_43^AssignVariableOp_44^AssignVariableOp_45^AssignVariableOp_46^AssignVariableOp_47^AssignVariableOp_48^AssignVariableOp_49^AssignVariableOp_5^AssignVariableOp_50^AssignVariableOp_51^AssignVariableOp_52^AssignVariableOp_53^AssignVariableOp_54^AssignVariableOp_6^AssignVariableOp_7^AssignVariableOp_8^AssignVariableOp_9^NoOp"/device:CPU:0*
T0*
_output_shapes
: 2
Identity_55

Identity_56IdentityIdentity_55:output:0^AssignVariableOp^AssignVariableOp_1^AssignVariableOp_10^AssignVariableOp_11^AssignVariableOp_12^AssignVariableOp_13^AssignVariableOp_14^AssignVariableOp_15^AssignVariableOp_16^AssignVariableOp_17^AssignVariableOp_18^AssignVariableOp_19^AssignVariableOp_2^AssignVariableOp_20^AssignVariableOp_21^AssignVariableOp_22^AssignVariableOp_23^AssignVariableOp_24^AssignVariableOp_25^AssignVariableOp_26^AssignVariableOp_27^AssignVariableOp_28^AssignVariableOp_29^AssignVariableOp_3^AssignVariableOp_30^AssignVariableOp_31^AssignVariableOp_32^AssignVariableOp_33^AssignVariableOp_34^AssignVariableOp_35^AssignVariableOp_36^AssignVariableOp_37^AssignVariableOp_38^AssignVariableOp_39^AssignVariableOp_4^AssignVariableOp_40^AssignVariableOp_41^AssignVariableOp_42^AssignVariableOp_43^AssignVariableOp_44^AssignVariableOp_45^AssignVariableOp_46^AssignVariableOp_47^AssignVariableOp_48^AssignVariableOp_49^AssignVariableOp_5^AssignVariableOp_50^AssignVariableOp_51^AssignVariableOp_52^AssignVariableOp_53^AssignVariableOp_54^AssignVariableOp_6^AssignVariableOp_7^AssignVariableOp_8^AssignVariableOp_9*
T0*
_output_shapes
: 2
Identity_56"#
identity_56Identity_56:output:0*ѓ
_input_shapesс
о: :::::::::::::::::::::::::::::::::::::::::::::::::::::::2$
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
AssignVariableOp_54AssignVariableOp_542(
AssignVariableOp_6AssignVariableOp_62(
AssignVariableOp_7AssignVariableOp_72(
AssignVariableOp_8AssignVariableOp_82(
AssignVariableOp_9AssignVariableOp_9:C ?

_output_shapes
: 
%
_user_specified_namefile_prefix
ѕ	
ф
K__inference_discriminator_l3_layer_call_and_return_conditional_losses_71021

inputs"
matmul_readvariableop_resource#
biasadd_readvariableop_resource
identityЂBiasAdd/ReadVariableOpЂMatMul/ReadVariableOp
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:*
dtype02
MatMul/ReadVariableOps
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2
MatMul
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:*
dtype02
BiasAdd/ReadVariableOp
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2	
BiasAddX
ReluReluBiasAdd:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
Relu
IdentityIdentityRelu:activations:0^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*
T0*'
_output_shapes
:џџџџџџџџџ2

Identity"
identityIdentity:output:0*.
_input_shapes
:џџџџџџџџџ::20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:O K
'
_output_shapes
:џџџџџџџџџ
 
_user_specified_nameinputs
ѕ	
ф
K__inference_generator_output_layer_call_and_return_conditional_losses_70820

inputs"
matmul_readvariableop_resource#
biasadd_readvariableop_resource
identityЂBiasAdd/ReadVariableOpЂMatMul/ReadVariableOp
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:*
dtype02
MatMul/ReadVariableOps
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2
MatMul
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:*
dtype02
BiasAdd/ReadVariableOp
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2	
BiasAddX
ReluReluBiasAdd:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
Relu
IdentityIdentityRelu:activations:0^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*
T0*'
_output_shapes
:џџџџџџџџџ2

Identity"
identityIdentity:output:0*.
_input_shapes
:џџџџџџџџџ::20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:O K
'
_output_shapes
:џџџџџџџџџ
 
_user_specified_nameinputs
ѕ	
ф
K__inference_discriminator_l3_layer_call_and_return_conditional_losses_71936

inputs"
matmul_readvariableop_resource#
biasadd_readvariableop_resource
identityЂBiasAdd/ReadVariableOpЂMatMul/ReadVariableOp
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:*
dtype02
MatMul/ReadVariableOps
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2
MatMul
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:*
dtype02
BiasAdd/ReadVariableOp
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2	
BiasAddX
ReluReluBiasAdd:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
Relu
IdentityIdentityRelu:activations:0^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*
T0*'
_output_shapes
:џџџџџџџџџ2

Identity"
identityIdentity:output:0*.
_input_shapes
:џџџџџџџџџ::20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:O K
'
_output_shapes
:џџџџџџџџџ
 
_user_specified_nameinputs
ы

0__inference_generator_output_layer_call_fn_71885

inputs
unknown
	unknown_0
identityЂStatefulPartitionedCallћ
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:џџџџџџџџџ*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8 *T
fORM
K__inference_generator_output_layer_call_and_return_conditional_losses_708202
StatefulPartitionedCall
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*'
_output_shapes
:џџџџџџџџџ2

Identity"
identityIdentity:output:0*.
_input_shapes
:џџџџџџџџџ::22
StatefulPartitionedCallStatefulPartitionedCall:O K
'
_output_shapes
:џџџџџџџџџ
 
_user_specified_nameinputs
ы

0__inference_discriminator_l3_layer_call_fn_71945

inputs
unknown
	unknown_0
identityЂStatefulPartitionedCallћ
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:џџџџџџџџџ*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8 *T
fORM
K__inference_discriminator_l3_layer_call_and_return_conditional_losses_710212
StatefulPartitionedCall
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*'
_output_shapes
:џџџџџџџџџ2

Identity"
identityIdentity:output:0*.
_input_shapes
:џџџџџџџџџ::22
StatefulPartitionedCallStatefulPartitionedCall:O K
'
_output_shapes
:џџџџџџџџџ
 
_user_specified_nameinputs


H__inference_discriminator_layer_call_and_return_conditional_losses_71065
discriminator_input
discriminator_l1_70978
discriminator_l1_70980
discriminator_l2_71005
discriminator_l2_71007
discriminator_l3_71032
discriminator_l3_71034
discriminator_output_71059
discriminator_output_71061
identityЂ(discriminator_l1/StatefulPartitionedCallЂ(discriminator_l2/StatefulPartitionedCallЂ(discriminator_l3/StatefulPartitionedCallЂ,discriminator_output/StatefulPartitionedCallЦ
(discriminator_l1/StatefulPartitionedCallStatefulPartitionedCalldiscriminator_inputdiscriminator_l1_70978discriminator_l1_70980*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:џџџџџџџџџ
*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8 *T
fORM
K__inference_discriminator_l1_layer_call_and_return_conditional_losses_709672*
(discriminator_l1/StatefulPartitionedCallф
(discriminator_l2/StatefulPartitionedCallStatefulPartitionedCall1discriminator_l1/StatefulPartitionedCall:output:0discriminator_l2_71005discriminator_l2_71007*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:џџџџџџџџџ*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8 *T
fORM
K__inference_discriminator_l2_layer_call_and_return_conditional_losses_709942*
(discriminator_l2/StatefulPartitionedCallф
(discriminator_l3/StatefulPartitionedCallStatefulPartitionedCall1discriminator_l2/StatefulPartitionedCall:output:0discriminator_l3_71032discriminator_l3_71034*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:џџџџџџџџџ*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8 *T
fORM
K__inference_discriminator_l3_layer_call_and_return_conditional_losses_710212*
(discriminator_l3/StatefulPartitionedCallј
,discriminator_output/StatefulPartitionedCallStatefulPartitionedCall1discriminator_l3/StatefulPartitionedCall:output:0discriminator_output_71059discriminator_output_71061*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:џџџџџџџџџ*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8 *X
fSRQ
O__inference_discriminator_output_layer_call_and_return_conditional_losses_710482.
,discriminator_output/StatefulPartitionedCallЙ
IdentityIdentity5discriminator_output/StatefulPartitionedCall:output:0)^discriminator_l1/StatefulPartitionedCall)^discriminator_l2/StatefulPartitionedCall)^discriminator_l3/StatefulPartitionedCall-^discriminator_output/StatefulPartitionedCall*
T0*'
_output_shapes
:џџџџџџџџџ2

Identity"
identityIdentity:output:0*F
_input_shapes5
3:џџџџџџџџџ::::::::2T
(discriminator_l1/StatefulPartitionedCall(discriminator_l1/StatefulPartitionedCall2T
(discriminator_l2/StatefulPartitionedCall(discriminator_l2/StatefulPartitionedCall2T
(discriminator_l3/StatefulPartitionedCall(discriminator_l3/StatefulPartitionedCall2\
,discriminator_output/StatefulPartitionedCall,discriminator_output/StatefulPartitionedCall:\ X
'
_output_shapes
:џџџџџџџџџ
-
_user_specified_namediscriminator_input
Р

б
B__inference_cgan_31_layer_call_and_return_conditional_losses_71331

inputs
discriminator_71313
discriminator_71315
discriminator_71317
discriminator_71319
discriminator_71321
discriminator_71323
discriminator_71325
discriminator_71327
identityЂ%discriminator/StatefulPartitionedCallД
%discriminator/StatefulPartitionedCallStatefulPartitionedCallinputsdiscriminator_71313discriminator_71315discriminator_71317discriminator_71319discriminator_71321discriminator_71323discriminator_71325discriminator_71327*
Tin
2	*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:џџџџџџџџџ**
_read_only_resource_inputs

*-
config_proto

CPU

GPU 2J 8 *Q
fLRJ
H__inference_discriminator_layer_call_and_return_conditional_losses_711612'
%discriminator/StatefulPartitionedCallЊ
IdentityIdentity.discriminator/StatefulPartitionedCall:output:0&^discriminator/StatefulPartitionedCall*
T0*'
_output_shapes
:џџџџџџџџџ2

Identity"
identityIdentity:output:0*F
_input_shapes5
3:џџџџџџџџџ::::::::2N
%discriminator/StatefulPartitionedCall%discriminator/StatefulPartitionedCall:O K
'
_output_shapes
:џџџџџџџџџ
 
_user_specified_nameinputs


H__inference_discriminator_layer_call_and_return_conditional_losses_71089
discriminator_input
discriminator_l1_71068
discriminator_l1_71070
discriminator_l2_71073
discriminator_l2_71075
discriminator_l3_71078
discriminator_l3_71080
discriminator_output_71083
discriminator_output_71085
identityЂ(discriminator_l1/StatefulPartitionedCallЂ(discriminator_l2/StatefulPartitionedCallЂ(discriminator_l3/StatefulPartitionedCallЂ,discriminator_output/StatefulPartitionedCallЦ
(discriminator_l1/StatefulPartitionedCallStatefulPartitionedCalldiscriminator_inputdiscriminator_l1_71068discriminator_l1_71070*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:џџџџџџџџџ
*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8 *T
fORM
K__inference_discriminator_l1_layer_call_and_return_conditional_losses_709672*
(discriminator_l1/StatefulPartitionedCallф
(discriminator_l2/StatefulPartitionedCallStatefulPartitionedCall1discriminator_l1/StatefulPartitionedCall:output:0discriminator_l2_71073discriminator_l2_71075*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:џџџџџџџџџ*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8 *T
fORM
K__inference_discriminator_l2_layer_call_and_return_conditional_losses_709942*
(discriminator_l2/StatefulPartitionedCallф
(discriminator_l3/StatefulPartitionedCallStatefulPartitionedCall1discriminator_l2/StatefulPartitionedCall:output:0discriminator_l3_71078discriminator_l3_71080*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:џџџџџџџџџ*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8 *T
fORM
K__inference_discriminator_l3_layer_call_and_return_conditional_losses_710212*
(discriminator_l3/StatefulPartitionedCallј
,discriminator_output/StatefulPartitionedCallStatefulPartitionedCall1discriminator_l3/StatefulPartitionedCall:output:0discriminator_output_71083discriminator_output_71085*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:џџџџџџџџџ*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8 *X
fSRQ
O__inference_discriminator_output_layer_call_and_return_conditional_losses_710482.
,discriminator_output/StatefulPartitionedCallЙ
IdentityIdentity5discriminator_output/StatefulPartitionedCall:output:0)^discriminator_l1/StatefulPartitionedCall)^discriminator_l2/StatefulPartitionedCall)^discriminator_l3/StatefulPartitionedCall-^discriminator_output/StatefulPartitionedCall*
T0*'
_output_shapes
:џџџџџџџџџ2

Identity"
identityIdentity:output:0*F
_input_shapes5
3:џџџџџџџџџ::::::::2T
(discriminator_l1/StatefulPartitionedCall(discriminator_l1/StatefulPartitionedCall2T
(discriminator_l2/StatefulPartitionedCall(discriminator_l2/StatefulPartitionedCall2T
(discriminator_l3/StatefulPartitionedCall(discriminator_l3/StatefulPartitionedCall2\
,discriminator_output/StatefulPartitionedCall,discriminator_output/StatefulPartitionedCall:\ X
'
_output_shapes
:џџџџџџџџџ
-
_user_specified_namediscriminator_input

ж
'__inference_cgan_31_layer_call_fn_71593

inputs
unknown
	unknown_0
	unknown_1
	unknown_2
	unknown_3
	unknown_4
	unknown_5
	unknown_6
identityЂStatefulPartitionedCallР
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6*
Tin
2	*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:џџџџџџџџџ**
_read_only_resource_inputs

*-
config_proto

CPU

GPU 2J 8 *K
fFRD
B__inference_cgan_31_layer_call_and_return_conditional_losses_713312
StatefulPartitionedCall
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*'
_output_shapes
:џџџџџџџџџ2

Identity"
identityIdentity:output:0*F
_input_shapes5
3:џџџџџџџџџ::::::::22
StatefulPartitionedCallStatefulPartitionedCall:O K
'
_output_shapes
:џџџџџџџџџ
 
_user_specified_nameinputs
ё	
р
G__inference_generator_l1_layer_call_and_return_conditional_losses_70739

inputs"
matmul_readvariableop_resource#
biasadd_readvariableop_resource
identityЂBiasAdd/ReadVariableOpЂMatMul/ReadVariableOp
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:
*
dtype02
MatMul/ReadVariableOps
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ
2
MatMul
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:
*
dtype02
BiasAdd/ReadVariableOp
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ
2	
BiasAddX
ReluReluBiasAdd:output:0*
T0*'
_output_shapes
:џџџџџџџџџ
2
Relu
IdentityIdentityRelu:activations:0^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*
T0*'
_output_shapes
:џџџџџџџџџ
2

Identity"
identityIdentity:output:0*.
_input_shapes
:џџџџџџџџџ::20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:O K
'
_output_shapes
:џџџџџџџџџ
 
_user_specified_nameinputs
Ѕ
м
-__inference_discriminator_layer_call_fn_71784

inputs
unknown
	unknown_0
	unknown_1
	unknown_2
	unknown_3
	unknown_4
	unknown_5
	unknown_6
identityЂStatefulPartitionedCallЦ
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6*
Tin
2	*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:џџџџџџџџџ**
_read_only_resource_inputs

*-
config_proto

CPU

GPU 2J 8 *Q
fLRJ
H__inference_discriminator_layer_call_and_return_conditional_losses_711162
StatefulPartitionedCall
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*'
_output_shapes
:џџџџџџџџџ2

Identity"
identityIdentity:output:0*F
_input_shapes5
3:џџџџџџџџџ::::::::22
StatefulPartitionedCallStatefulPartitionedCall:O K
'
_output_shapes
:џџџџџџџџџ
 
_user_specified_nameinputs"БL
saver_filename:0StatefulPartitionedCall_1:0StatefulPartitionedCall_28"
saved_model_main_op

NoOp*>
__saved_model_init_op%#
__saved_model_init_op

NoOp*Ћ
serving_default
;
input_10
serving_default_input_1:0џџџџџџџџџ<
output_10
StatefulPartitionedCall:0џџџџџџџџџtensorflow/serving/predict:др
с
	optimizer
	generator
discriminator
generator_optimizer
discriminator_optimizer
loss
regularization_losses
trainable_variables
	variables
	keras_api
	
signatures
­_default_save_signature
Ў__call__
+Џ&call_and_return_all_conditional_losses"Б
_tf_keras_model{"class_name": "CGAN", "name": "cgan_31", "trainable": true, "expects_training_arg": true, "dtype": "float32", "batch_input_shape": null, "must_restore_from_config": false, "config": {"layer was saved without config": true}, "is_graph_network": false, "keras_version": "2.4.0", "backend": "tensorflow", "model_config": {"class_name": "CGAN"}, "training_config": {"loss": null, "metrics": null, "weighted_metrics": null, "loss_weights": null, "optimizer_config": {"class_name": "Adam", "config": {"name": "Adam", "learning_rate": 0.009999999776482582, "decay": 0.0, "beta_1": 0.8999999761581421, "beta_2": 0.9990000128746033, "epsilon": 1e-07, "amsgrad": false}}}}


iter

beta_1

beta_2
	decay
learning_rate!m"m#m$m%m&m'm(m)m*m+m,m-m.m/m0m!v"v#v$v %vЁ&vЂ'vЃ(vЄ)vЅ*vІ+vЇ,vЈ-vЉ.vЊ/vЋ0vЌ"
	optimizer
Р,
layer-0
layer_with_weights-0
layer-1
layer_with_weights-1
layer-2
layer_with_weights-2
layer-3
layer_with_weights-3
layer-4
regularization_losses
trainable_variables
	variables
	keras_api
А__call__
+Б&call_and_return_all_conditional_losses"*
_tf_keras_networkъ){"class_name": "Functional", "name": "generator", "trainable": true, "expects_training_arg": true, "dtype": "float32", "batch_input_shape": null, "must_restore_from_config": false, "config": {"name": "generator", "layers": [{"class_name": "InputLayer", "config": {"batch_input_shape": {"class_name": "__tuple__", "items": [null, 4]}, "dtype": "float32", "sparse": false, "ragged": false, "name": "generator_input"}, "name": "generator_input", "inbound_nodes": []}, {"class_name": "Dense", "config": {"name": "generator_l1", "trainable": true, "dtype": "float32", "units": 10, "activation": "relu", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "name": "generator_l1", "inbound_nodes": [[["generator_input", 0, 0, {}]]]}, {"class_name": "Dense", "config": {"name": "generator_l2", "trainable": true, "dtype": "float32", "units": 8, "activation": "relu", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "name": "generator_l2", "inbound_nodes": [[["generator_l1", 0, 0, {}]]]}, {"class_name": "Dense", "config": {"name": "generator_l3", "trainable": true, "dtype": "float32", "units": 6, "activation": "relu", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "name": "generator_l3", "inbound_nodes": [[["generator_l2", 0, 0, {}]]]}, {"class_name": "Dense", "config": {"name": "generator_output", "trainable": true, "dtype": "float32", "units": 12, "activation": "relu", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "name": "generator_output", "inbound_nodes": [[["generator_l3", 0, 0, {}]]]}], "input_layers": [["generator_input", 0, 0]], "output_layers": [["generator_output", 0, 0]]}, "input_spec": [{"class_name": "InputSpec", "config": {"dtype": null, "shape": {"class_name": "__tuple__", "items": [null, 4]}, "ndim": 2, "max_ndim": null, "min_ndim": null, "axes": {}}}], "build_input_shape": {"class_name": "TensorShape", "items": [null, 4]}, "is_graph_network": true, "keras_version": "2.4.0", "backend": "tensorflow", "model_config": {"class_name": "Functional", "config": {"name": "generator", "layers": [{"class_name": "InputLayer", "config": {"batch_input_shape": {"class_name": "__tuple__", "items": [null, 4]}, "dtype": "float32", "sparse": false, "ragged": false, "name": "generator_input"}, "name": "generator_input", "inbound_nodes": []}, {"class_name": "Dense", "config": {"name": "generator_l1", "trainable": true, "dtype": "float32", "units": 10, "activation": "relu", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "name": "generator_l1", "inbound_nodes": [[["generator_input", 0, 0, {}]]]}, {"class_name": "Dense", "config": {"name": "generator_l2", "trainable": true, "dtype": "float32", "units": 8, "activation": "relu", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "name": "generator_l2", "inbound_nodes": [[["generator_l1", 0, 0, {}]]]}, {"class_name": "Dense", "config": {"name": "generator_l3", "trainable": true, "dtype": "float32", "units": 6, "activation": "relu", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "name": "generator_l3", "inbound_nodes": [[["generator_l2", 0, 0, {}]]]}, {"class_name": "Dense", "config": {"name": "generator_output", "trainable": true, "dtype": "float32", "units": 12, "activation": "relu", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "name": "generator_output", "inbound_nodes": [[["generator_l3", 0, 0, {}]]]}], "input_layers": [["generator_input", 0, 0]], "output_layers": [["generator_output", 0, 0]]}}}
д-
layer-0
layer_with_weights-0
layer-1
layer_with_weights-1
layer-2
layer_with_weights-2
layer-3
layer_with_weights-3
layer-4
regularization_losses
trainable_variables
	variables
 	keras_api
В__call__
+Г&call_and_return_all_conditional_losses"+
_tf_keras_networkў*{"class_name": "Functional", "name": "discriminator", "trainable": true, "expects_training_arg": true, "dtype": "float32", "batch_input_shape": null, "must_restore_from_config": false, "config": {"name": "discriminator", "layers": [{"class_name": "InputLayer", "config": {"batch_input_shape": {"class_name": "__tuple__", "items": [null, 14]}, "dtype": "float32", "sparse": false, "ragged": false, "name": "discriminator_input"}, "name": "discriminator_input", "inbound_nodes": []}, {"class_name": "Dense", "config": {"name": "discriminator_l1", "trainable": true, "dtype": "float32", "units": 10, "activation": "relu", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "name": "discriminator_l1", "inbound_nodes": [[["discriminator_input", 0, 0, {}]]]}, {"class_name": "Dense", "config": {"name": "discriminator_l2", "trainable": true, "dtype": "float32", "units": 8, "activation": "relu", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "name": "discriminator_l2", "inbound_nodes": [[["discriminator_l1", 0, 0, {}]]]}, {"class_name": "Dense", "config": {"name": "discriminator_l3", "trainable": true, "dtype": "float32", "units": 6, "activation": "relu", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "name": "discriminator_l3", "inbound_nodes": [[["discriminator_l2", 0, 0, {}]]]}, {"class_name": "Dense", "config": {"name": "discriminator_output", "trainable": true, "dtype": "float32", "units": 1, "activation": "sigmoid", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "name": "discriminator_output", "inbound_nodes": [[["discriminator_l3", 0, 0, {}]]]}], "input_layers": [["discriminator_input", 0, 0]], "output_layers": [["discriminator_output", 0, 0]]}, "input_spec": [{"class_name": "InputSpec", "config": {"dtype": null, "shape": {"class_name": "__tuple__", "items": [null, 14]}, "ndim": 2, "max_ndim": null, "min_ndim": null, "axes": {}}}], "build_input_shape": {"class_name": "TensorShape", "items": [null, 14]}, "is_graph_network": true, "keras_version": "2.4.0", "backend": "tensorflow", "model_config": {"class_name": "Functional", "config": {"name": "discriminator", "layers": [{"class_name": "InputLayer", "config": {"batch_input_shape": {"class_name": "__tuple__", "items": [null, 14]}, "dtype": "float32", "sparse": false, "ragged": false, "name": "discriminator_input"}, "name": "discriminator_input", "inbound_nodes": []}, {"class_name": "Dense", "config": {"name": "discriminator_l1", "trainable": true, "dtype": "float32", "units": 10, "activation": "relu", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "name": "discriminator_l1", "inbound_nodes": [[["discriminator_input", 0, 0, {}]]]}, {"class_name": "Dense", "config": {"name": "discriminator_l2", "trainable": true, "dtype": "float32", "units": 8, "activation": "relu", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "name": "discriminator_l2", "inbound_nodes": [[["discriminator_l1", 0, 0, {}]]]}, {"class_name": "Dense", "config": {"name": "discriminator_l3", "trainable": true, "dtype": "float32", "units": 6, "activation": "relu", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "name": "discriminator_l3", "inbound_nodes": [[["discriminator_l2", 0, 0, {}]]]}, {"class_name": "Dense", "config": {"name": "discriminator_output", "trainable": true, "dtype": "float32", "units": 1, "activation": "sigmoid", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "name": "discriminator_output", "inbound_nodes": [[["discriminator_l3", 0, 0, {}]]]}], "input_layers": [["discriminator_input", 0, 0]], "output_layers": [["discriminator_output", 0, 0]]}}}
 "
trackable_dict_wrapper
 "
trackable_list_wrapper

!0
"1
#2
$3
%4
&5
'6
(7
)8
*9
+10
,11
-12
.13
/14
015"
trackable_list_wrapper

!0
"1
#2
$3
%4
&5
'6
(7
)8
*9
+10
,11
-12
.13
/14
015"
trackable_list_wrapper
Ю
1layer_regularization_losses

2layers
3layer_metrics
4non_trainable_variables
5metrics
regularization_losses
trainable_variables
	variables
Ў__call__
­_default_save_signature
+Џ&call_and_return_all_conditional_losses
'Џ"call_and_return_conditional_losses"
_generic_user_object
-
Дserving_default"
signature_map
:	 (2	Adam/iter
: (2Adam/beta_1
: (2Adam/beta_2
: (2
Adam/decay
: (2Adam/learning_rate
љ"і
_tf_keras_input_layerж{"class_name": "InputLayer", "name": "generator_input", "dtype": "float32", "sparse": false, "ragged": false, "batch_input_shape": {"class_name": "__tuple__", "items": [null, 4]}, "config": {"batch_input_shape": {"class_name": "__tuple__", "items": [null, 4]}, "dtype": "float32", "sparse": false, "ragged": false, "name": "generator_input"}}
њ

!kernel
"bias
6regularization_losses
7trainable_variables
8	variables
9	keras_api
Е__call__
+Ж&call_and_return_all_conditional_losses"г
_tf_keras_layerЙ{"class_name": "Dense", "name": "generator_l1", "trainable": true, "expects_training_arg": false, "dtype": "float32", "batch_input_shape": null, "stateful": false, "must_restore_from_config": false, "config": {"name": "generator_l1", "trainable": true, "dtype": "float32", "units": 10, "activation": "relu", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "input_spec": {"class_name": "InputSpec", "config": {"dtype": null, "shape": null, "ndim": null, "max_ndim": null, "min_ndim": 2, "axes": {"-1": 4}}}, "build_input_shape": {"class_name": "TensorShape", "items": [null, 4]}}
ћ

#kernel
$bias
:regularization_losses
;trainable_variables
<	variables
=	keras_api
З__call__
+И&call_and_return_all_conditional_losses"д
_tf_keras_layerК{"class_name": "Dense", "name": "generator_l2", "trainable": true, "expects_training_arg": false, "dtype": "float32", "batch_input_shape": null, "stateful": false, "must_restore_from_config": false, "config": {"name": "generator_l2", "trainable": true, "dtype": "float32", "units": 8, "activation": "relu", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "input_spec": {"class_name": "InputSpec", "config": {"dtype": null, "shape": null, "ndim": null, "max_ndim": null, "min_ndim": 2, "axes": {"-1": 10}}}, "build_input_shape": {"class_name": "TensorShape", "items": [null, 10]}}
љ

%kernel
&bias
>regularization_losses
?trainable_variables
@	variables
A	keras_api
Й__call__
+К&call_and_return_all_conditional_losses"в
_tf_keras_layerИ{"class_name": "Dense", "name": "generator_l3", "trainable": true, "expects_training_arg": false, "dtype": "float32", "batch_input_shape": null, "stateful": false, "must_restore_from_config": false, "config": {"name": "generator_l3", "trainable": true, "dtype": "float32", "units": 6, "activation": "relu", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "input_spec": {"class_name": "InputSpec", "config": {"dtype": null, "shape": null, "ndim": null, "max_ndim": null, "min_ndim": 2, "axes": {"-1": 8}}}, "build_input_shape": {"class_name": "TensorShape", "items": [null, 8]}}


'kernel
(bias
Bregularization_losses
Ctrainable_variables
D	variables
E	keras_api
Л__call__
+М&call_and_return_all_conditional_losses"л
_tf_keras_layerС{"class_name": "Dense", "name": "generator_output", "trainable": true, "expects_training_arg": false, "dtype": "float32", "batch_input_shape": null, "stateful": false, "must_restore_from_config": false, "config": {"name": "generator_output", "trainable": true, "dtype": "float32", "units": 12, "activation": "relu", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "input_spec": {"class_name": "InputSpec", "config": {"dtype": null, "shape": null, "ndim": null, "max_ndim": null, "min_ndim": 2, "axes": {"-1": 6}}}, "build_input_shape": {"class_name": "TensorShape", "items": [null, 6]}}
 "
trackable_list_wrapper
X
!0
"1
#2
$3
%4
&5
'6
(7"
trackable_list_wrapper
X
!0
"1
#2
$3
%4
&5
'6
(7"
trackable_list_wrapper
А
Flayer_regularization_losses

Glayers
Hlayer_metrics
Inon_trainable_variables
Jmetrics
regularization_losses
trainable_variables
	variables
А__call__
+Б&call_and_return_all_conditional_losses
'Б"call_and_return_conditional_losses"
_generic_user_object
"
_tf_keras_input_layerр{"class_name": "InputLayer", "name": "discriminator_input", "dtype": "float32", "sparse": false, "ragged": false, "batch_input_shape": {"class_name": "__tuple__", "items": [null, 14]}, "config": {"batch_input_shape": {"class_name": "__tuple__", "items": [null, 14]}, "dtype": "float32", "sparse": false, "ragged": false, "name": "discriminator_input"}}


)kernel
*bias
Kregularization_losses
Ltrainable_variables
M	variables
N	keras_api
Н__call__
+О&call_and_return_all_conditional_losses"н
_tf_keras_layerУ{"class_name": "Dense", "name": "discriminator_l1", "trainable": true, "expects_training_arg": false, "dtype": "float32", "batch_input_shape": null, "stateful": false, "must_restore_from_config": false, "config": {"name": "discriminator_l1", "trainable": true, "dtype": "float32", "units": 10, "activation": "relu", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "input_spec": {"class_name": "InputSpec", "config": {"dtype": null, "shape": null, "ndim": null, "max_ndim": null, "min_ndim": 2, "axes": {"-1": 14}}}, "build_input_shape": {"class_name": "TensorShape", "items": [null, 14]}}


+kernel
,bias
Oregularization_losses
Ptrainable_variables
Q	variables
R	keras_api
П__call__
+Р&call_and_return_all_conditional_losses"м
_tf_keras_layerТ{"class_name": "Dense", "name": "discriminator_l2", "trainable": true, "expects_training_arg": false, "dtype": "float32", "batch_input_shape": null, "stateful": false, "must_restore_from_config": false, "config": {"name": "discriminator_l2", "trainable": true, "dtype": "float32", "units": 8, "activation": "relu", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "input_spec": {"class_name": "InputSpec", "config": {"dtype": null, "shape": null, "ndim": null, "max_ndim": null, "min_ndim": 2, "axes": {"-1": 10}}}, "build_input_shape": {"class_name": "TensorShape", "items": [null, 10]}}


-kernel
.bias
Sregularization_losses
Ttrainable_variables
U	variables
V	keras_api
С__call__
+Т&call_and_return_all_conditional_losses"к
_tf_keras_layerР{"class_name": "Dense", "name": "discriminator_l3", "trainable": true, "expects_training_arg": false, "dtype": "float32", "batch_input_shape": null, "stateful": false, "must_restore_from_config": false, "config": {"name": "discriminator_l3", "trainable": true, "dtype": "float32", "units": 6, "activation": "relu", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "input_spec": {"class_name": "InputSpec", "config": {"dtype": null, "shape": null, "ndim": null, "max_ndim": null, "min_ndim": 2, "axes": {"-1": 8}}}, "build_input_shape": {"class_name": "TensorShape", "items": [null, 8]}}


/kernel
0bias
Wregularization_losses
Xtrainable_variables
Y	variables
Z	keras_api
У__call__
+Ф&call_and_return_all_conditional_losses"х
_tf_keras_layerЫ{"class_name": "Dense", "name": "discriminator_output", "trainable": true, "expects_training_arg": false, "dtype": "float32", "batch_input_shape": null, "stateful": false, "must_restore_from_config": false, "config": {"name": "discriminator_output", "trainable": true, "dtype": "float32", "units": 1, "activation": "sigmoid", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "input_spec": {"class_name": "InputSpec", "config": {"dtype": null, "shape": null, "ndim": null, "max_ndim": null, "min_ndim": 2, "axes": {"-1": 6}}}, "build_input_shape": {"class_name": "TensorShape", "items": [null, 6]}}
 "
trackable_list_wrapper
X
)0
*1
+2
,3
-4
.5
/6
07"
trackable_list_wrapper
X
)0
*1
+2
,3
-4
.5
/6
07"
trackable_list_wrapper
А
[layer_regularization_losses

\layers
]layer_metrics
^non_trainable_variables
_metrics
regularization_losses
trainable_variables
	variables
В__call__
+Г&call_and_return_all_conditional_losses
'Г"call_and_return_conditional_losses"
_generic_user_object
%:#
2generator_l1/kernel
:
2generator_l1/bias
%:#
2generator_l2/kernel
:2generator_l2/bias
%:#2generator_l3/kernel
:2generator_l3/bias
):'2generator_output/kernel
#:!2generator_output/bias
):'
2discriminator_l1/kernel
#:!
2discriminator_l1/bias
):'
2discriminator_l2/kernel
#:!2discriminator_l2/bias
):'2discriminator_l3/kernel
#:!2discriminator_l3/bias
-:+2discriminator_output/kernel
':%2discriminator_output/bias
 "
trackable_list_wrapper
.
0
1"
trackable_list_wrapper
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
'
`0"
trackable_list_wrapper
 "
trackable_list_wrapper
.
!0
"1"
trackable_list_wrapper
.
!0
"1"
trackable_list_wrapper
А

alayers
blayer_regularization_losses
clayer_metrics
dnon_trainable_variables
emetrics
6regularization_losses
7trainable_variables
8	variables
Е__call__
+Ж&call_and_return_all_conditional_losses
'Ж"call_and_return_conditional_losses"
_generic_user_object
 "
trackable_list_wrapper
.
#0
$1"
trackable_list_wrapper
.
#0
$1"
trackable_list_wrapper
А

flayers
glayer_regularization_losses
hlayer_metrics
inon_trainable_variables
jmetrics
:regularization_losses
;trainable_variables
<	variables
З__call__
+И&call_and_return_all_conditional_losses
'И"call_and_return_conditional_losses"
_generic_user_object
 "
trackable_list_wrapper
.
%0
&1"
trackable_list_wrapper
.
%0
&1"
trackable_list_wrapper
А

klayers
llayer_regularization_losses
mlayer_metrics
nnon_trainable_variables
ometrics
>regularization_losses
?trainable_variables
@	variables
Й__call__
+К&call_and_return_all_conditional_losses
'К"call_and_return_conditional_losses"
_generic_user_object
 "
trackable_list_wrapper
.
'0
(1"
trackable_list_wrapper
.
'0
(1"
trackable_list_wrapper
А

players
qlayer_regularization_losses
rlayer_metrics
snon_trainable_variables
tmetrics
Bregularization_losses
Ctrainable_variables
D	variables
Л__call__
+М&call_and_return_all_conditional_losses
'М"call_and_return_conditional_losses"
_generic_user_object
 "
trackable_list_wrapper
C
0
1
2
3
4"
trackable_list_wrapper
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
.
)0
*1"
trackable_list_wrapper
.
)0
*1"
trackable_list_wrapper
А

ulayers
vlayer_regularization_losses
wlayer_metrics
xnon_trainable_variables
ymetrics
Kregularization_losses
Ltrainable_variables
M	variables
Н__call__
+О&call_and_return_all_conditional_losses
'О"call_and_return_conditional_losses"
_generic_user_object
 "
trackable_list_wrapper
.
+0
,1"
trackable_list_wrapper
.
+0
,1"
trackable_list_wrapper
А

zlayers
{layer_regularization_losses
|layer_metrics
}non_trainable_variables
~metrics
Oregularization_losses
Ptrainable_variables
Q	variables
П__call__
+Р&call_and_return_all_conditional_losses
'Р"call_and_return_conditional_losses"
_generic_user_object
 "
trackable_list_wrapper
.
-0
.1"
trackable_list_wrapper
.
-0
.1"
trackable_list_wrapper
Д

layers
 layer_regularization_losses
layer_metrics
non_trainable_variables
metrics
Sregularization_losses
Ttrainable_variables
U	variables
С__call__
+Т&call_and_return_all_conditional_losses
'Т"call_and_return_conditional_losses"
_generic_user_object
 "
trackable_list_wrapper
.
/0
01"
trackable_list_wrapper
.
/0
01"
trackable_list_wrapper
Е
layers
 layer_regularization_losses
layer_metrics
non_trainable_variables
metrics
Wregularization_losses
Xtrainable_variables
Y	variables
У__call__
+Ф&call_and_return_all_conditional_losses
'Ф"call_and_return_conditional_losses"
_generic_user_object
 "
trackable_list_wrapper
C
0
1
2
3
4"
trackable_list_wrapper
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
П

total

count
	variables
	keras_api"
_tf_keras_metricj{"class_name": "Mean", "name": "loss", "dtype": "float32", "config": {"name": "loss", "dtype": "float32"}}
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
:  (2total
:  (2count
0
0
1"
trackable_list_wrapper
.
	variables"
_generic_user_object
*:(
2Adam/generator_l1/kernel/m
$:"
2Adam/generator_l1/bias/m
*:(
2Adam/generator_l2/kernel/m
$:"2Adam/generator_l2/bias/m
*:(2Adam/generator_l3/kernel/m
$:"2Adam/generator_l3/bias/m
.:,2Adam/generator_output/kernel/m
(:&2Adam/generator_output/bias/m
0:.
2 Adam_1/discriminator_l1/kernel/m
*:(
2Adam_1/discriminator_l1/bias/m
0:.
2 Adam_1/discriminator_l2/kernel/m
*:(2Adam_1/discriminator_l2/bias/m
0:.2 Adam_1/discriminator_l3/kernel/m
*:(2Adam_1/discriminator_l3/bias/m
4:22$Adam_1/discriminator_output/kernel/m
.:,2"Adam_1/discriminator_output/bias/m
*:(
2Adam/generator_l1/kernel/v
$:"
2Adam/generator_l1/bias/v
*:(
2Adam/generator_l2/kernel/v
$:"2Adam/generator_l2/bias/v
*:(2Adam/generator_l3/kernel/v
$:"2Adam/generator_l3/bias/v
.:,2Adam/generator_output/kernel/v
(:&2Adam/generator_output/bias/v
0:.
2 Adam_1/discriminator_l1/kernel/v
*:(
2Adam_1/discriminator_l1/bias/v
0:.
2 Adam_1/discriminator_l2/kernel/v
*:(2Adam_1/discriminator_l2/bias/v
0:.2 Adam_1/discriminator_l3/kernel/v
*:(2Adam_1/discriminator_l3/bias/v
4:22$Adam_1/discriminator_output/kernel/v
.:,2"Adam_1/discriminator_output/bias/v
о2л
 __inference__wrapped_model_70724Ж
В
FullArgSpec
args 
varargsjargs
varkw
 
defaults
 

kwonlyargs 
kwonlydefaults
 
annotationsЊ *&Ђ#
!
input_1џџџџџџџџџ
о2л
'__inference_cgan_31_layer_call_fn_71487
'__inference_cgan_31_layer_call_fn_71466
'__inference_cgan_31_layer_call_fn_71572
'__inference_cgan_31_layer_call_fn_71593Д
ЋВЇ
FullArgSpec)
args!
jself
jinputs

jtraining
varargs
 
varkw
 
defaults
p 

kwonlyargs 
kwonlydefaultsЊ 
annotationsЊ *
 
Ъ2Ч
B__inference_cgan_31_layer_call_and_return_conditional_losses_71445
B__inference_cgan_31_layer_call_and_return_conditional_losses_71519
B__inference_cgan_31_layer_call_and_return_conditional_losses_71413
B__inference_cgan_31_layer_call_and_return_conditional_losses_71551Д
ЋВЇ
FullArgSpec)
args!
jself
jinputs

jtraining
varargs
 
varkw
 
defaults
p 

kwonlyargs 
kwonlydefaultsЊ 
annotationsЊ *
 
ђ2я
)__inference_generator_layer_call_fn_71699
)__inference_generator_layer_call_fn_70952
)__inference_generator_layer_call_fn_70907
)__inference_generator_layer_call_fn_71678Р
ЗВГ
FullArgSpec1
args)&
jself
jinputs

jtraining
jmask
varargs
 
varkw
 
defaults
p 

 

kwonlyargs 
kwonlydefaultsЊ 
annotationsЊ *
 
о2л
D__inference_generator_layer_call_and_return_conditional_losses_71625
D__inference_generator_layer_call_and_return_conditional_losses_71657
D__inference_generator_layer_call_and_return_conditional_losses_70837
D__inference_generator_layer_call_and_return_conditional_losses_70861Р
ЗВГ
FullArgSpec1
args)&
jself
jinputs

jtraining
jmask
varargs
 
varkw
 
defaults
p 

 

kwonlyargs 
kwonlydefaultsЊ 
annotationsЊ *
 
2џ
-__inference_discriminator_layer_call_fn_71135
-__inference_discriminator_layer_call_fn_71805
-__inference_discriminator_layer_call_fn_71180
-__inference_discriminator_layer_call_fn_71784Р
ЗВГ
FullArgSpec1
args)&
jself
jinputs

jtraining
jmask
varargs
 
varkw
 
defaults
p 

 

kwonlyargs 
kwonlydefaultsЊ 
annotationsЊ *
 
ю2ы
H__inference_discriminator_layer_call_and_return_conditional_losses_71731
H__inference_discriminator_layer_call_and_return_conditional_losses_71089
H__inference_discriminator_layer_call_and_return_conditional_losses_71763
H__inference_discriminator_layer_call_and_return_conditional_losses_71065Р
ЗВГ
FullArgSpec1
args)&
jself
jinputs

jtraining
jmask
varargs
 
varkw
 
defaults
p 

 

kwonlyargs 
kwonlydefaultsЊ 
annotationsЊ *
 
ЪBЧ
#__inference_signature_wrapper_71381input_1"
В
FullArgSpec
args 
varargs
 
varkwjkwargs
defaults
 

kwonlyargs 
kwonlydefaults
 
annotationsЊ *
 
ж2г
,__inference_generator_l1_layer_call_fn_71825Ђ
В
FullArgSpec
args
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs 
kwonlydefaults
 
annotationsЊ *
 
ё2ю
G__inference_generator_l1_layer_call_and_return_conditional_losses_71816Ђ
В
FullArgSpec
args
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs 
kwonlydefaults
 
annotationsЊ *
 
ж2г
,__inference_generator_l2_layer_call_fn_71845Ђ
В
FullArgSpec
args
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs 
kwonlydefaults
 
annotationsЊ *
 
ё2ю
G__inference_generator_l2_layer_call_and_return_conditional_losses_71836Ђ
В
FullArgSpec
args
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs 
kwonlydefaults
 
annotationsЊ *
 
ж2г
,__inference_generator_l3_layer_call_fn_71865Ђ
В
FullArgSpec
args
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs 
kwonlydefaults
 
annotationsЊ *
 
ё2ю
G__inference_generator_l3_layer_call_and_return_conditional_losses_71856Ђ
В
FullArgSpec
args
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs 
kwonlydefaults
 
annotationsЊ *
 
к2з
0__inference_generator_output_layer_call_fn_71885Ђ
В
FullArgSpec
args
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs 
kwonlydefaults
 
annotationsЊ *
 
ѕ2ђ
K__inference_generator_output_layer_call_and_return_conditional_losses_71876Ђ
В
FullArgSpec
args
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs 
kwonlydefaults
 
annotationsЊ *
 
к2з
0__inference_discriminator_l1_layer_call_fn_71905Ђ
В
FullArgSpec
args
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs 
kwonlydefaults
 
annotationsЊ *
 
ѕ2ђ
K__inference_discriminator_l1_layer_call_and_return_conditional_losses_71896Ђ
В
FullArgSpec
args
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs 
kwonlydefaults
 
annotationsЊ *
 
к2з
0__inference_discriminator_l2_layer_call_fn_71925Ђ
В
FullArgSpec
args
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs 
kwonlydefaults
 
annotationsЊ *
 
ѕ2ђ
K__inference_discriminator_l2_layer_call_and_return_conditional_losses_71916Ђ
В
FullArgSpec
args
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs 
kwonlydefaults
 
annotationsЊ *
 
к2з
0__inference_discriminator_l3_layer_call_fn_71945Ђ
В
FullArgSpec
args
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs 
kwonlydefaults
 
annotationsЊ *
 
ѕ2ђ
K__inference_discriminator_l3_layer_call_and_return_conditional_losses_71936Ђ
В
FullArgSpec
args
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs 
kwonlydefaults
 
annotationsЊ *
 
о2л
4__inference_discriminator_output_layer_call_fn_71965Ђ
В
FullArgSpec
args
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs 
kwonlydefaults
 
annotationsЊ *
 
љ2і
O__inference_discriminator_output_layer_call_and_return_conditional_losses_71956Ђ
В
FullArgSpec
args
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs 
kwonlydefaults
 
annotationsЊ *
 
 __inference__wrapped_model_70724q)*+,-./00Ђ-
&Ђ#
!
input_1џџџџџџџџџ
Њ "3Њ0
.
output_1"
output_1џџџџџџџџџ­
B__inference_cgan_31_layer_call_and_return_conditional_losses_71413g)*+,-./04Ђ1
*Ђ'
!
input_1џџџџџџџџџ
p
Њ "%Ђ"

0џџџџџџџџџ
 ­
B__inference_cgan_31_layer_call_and_return_conditional_losses_71445g)*+,-./04Ђ1
*Ђ'
!
input_1џџџџџџџџџ
p 
Њ "%Ђ"

0џџџџџџџџџ
 Ќ
B__inference_cgan_31_layer_call_and_return_conditional_losses_71519f)*+,-./03Ђ0
)Ђ&
 
inputsџџџџџџџџџ
p
Њ "%Ђ"

0џџџџџџџџџ
 Ќ
B__inference_cgan_31_layer_call_and_return_conditional_losses_71551f)*+,-./03Ђ0
)Ђ&
 
inputsџџџџџџџџџ
p 
Њ "%Ђ"

0џџџџџџџџџ
 
'__inference_cgan_31_layer_call_fn_71466Z)*+,-./04Ђ1
*Ђ'
!
input_1џџџџџџџџџ
p
Њ "џџџџџџџџџ
'__inference_cgan_31_layer_call_fn_71487Z)*+,-./04Ђ1
*Ђ'
!
input_1џџџџџџџџџ
p 
Њ "џџџџџџџџџ
'__inference_cgan_31_layer_call_fn_71572Y)*+,-./03Ђ0
)Ђ&
 
inputsџџџџџџџџџ
p
Њ "џџџџџџџџџ
'__inference_cgan_31_layer_call_fn_71593Y)*+,-./03Ђ0
)Ђ&
 
inputsџџџџџџџџџ
p 
Њ "џџџџџџџџџЋ
K__inference_discriminator_l1_layer_call_and_return_conditional_losses_71896\)*/Ђ,
%Ђ"
 
inputsџџџџџџџџџ
Њ "%Ђ"

0џџџџџџџџџ

 
0__inference_discriminator_l1_layer_call_fn_71905O)*/Ђ,
%Ђ"
 
inputsџџџџџџџџџ
Њ "џџџџџџџџџ
Ћ
K__inference_discriminator_l2_layer_call_and_return_conditional_losses_71916\+,/Ђ,
%Ђ"
 
inputsџџџџџџџџџ

Њ "%Ђ"

0џџџџџџџџџ
 
0__inference_discriminator_l2_layer_call_fn_71925O+,/Ђ,
%Ђ"
 
inputsџџџџџџџџџ

Њ "џџџџџџџџџЋ
K__inference_discriminator_l3_layer_call_and_return_conditional_losses_71936\-./Ђ,
%Ђ"
 
inputsџџџџџџџџџ
Њ "%Ђ"

0џџџџџџџџџ
 
0__inference_discriminator_l3_layer_call_fn_71945O-./Ђ,
%Ђ"
 
inputsџџџџџџџџџ
Њ "џџџџџџџџџУ
H__inference_discriminator_layer_call_and_return_conditional_losses_71065w)*+,-./0DЂA
:Ђ7
-*
discriminator_inputџџџџџџџџџ
p

 
Њ "%Ђ"

0џџџџџџџџџ
 У
H__inference_discriminator_layer_call_and_return_conditional_losses_71089w)*+,-./0DЂA
:Ђ7
-*
discriminator_inputџџџџџџџџџ
p 

 
Њ "%Ђ"

0џџџџџџџџџ
 Ж
H__inference_discriminator_layer_call_and_return_conditional_losses_71731j)*+,-./07Ђ4
-Ђ*
 
inputsџџџџџџџџџ
p

 
Њ "%Ђ"

0џџџџџџџџџ
 Ж
H__inference_discriminator_layer_call_and_return_conditional_losses_71763j)*+,-./07Ђ4
-Ђ*
 
inputsџџџџџџџџџ
p 

 
Њ "%Ђ"

0џџџџџџџџџ
 
-__inference_discriminator_layer_call_fn_71135j)*+,-./0DЂA
:Ђ7
-*
discriminator_inputџџџџџџџџџ
p

 
Њ "џџџџџџџџџ
-__inference_discriminator_layer_call_fn_71180j)*+,-./0DЂA
:Ђ7
-*
discriminator_inputџџџџџџџџџ
p 

 
Њ "џџџџџџџџџ
-__inference_discriminator_layer_call_fn_71784])*+,-./07Ђ4
-Ђ*
 
inputsџџџџџџџџџ
p

 
Њ "џџџџџџџџџ
-__inference_discriminator_layer_call_fn_71805])*+,-./07Ђ4
-Ђ*
 
inputsџџџџџџџџџ
p 

 
Њ "џџџџџџџџџЏ
O__inference_discriminator_output_layer_call_and_return_conditional_losses_71956\/0/Ђ,
%Ђ"
 
inputsџџџџџџџџџ
Њ "%Ђ"

0џџџџџџџџџ
 
4__inference_discriminator_output_layer_call_fn_71965O/0/Ђ,
%Ђ"
 
inputsџџџџџџџџџ
Њ "џџџџџџџџџЇ
G__inference_generator_l1_layer_call_and_return_conditional_losses_71816\!"/Ђ,
%Ђ"
 
inputsџџџџџџџџџ
Њ "%Ђ"

0џџџџџџџџџ

 
,__inference_generator_l1_layer_call_fn_71825O!"/Ђ,
%Ђ"
 
inputsџџџџџџџџџ
Њ "џџџџџџџџџ
Ї
G__inference_generator_l2_layer_call_and_return_conditional_losses_71836\#$/Ђ,
%Ђ"
 
inputsџџџџџџџџџ

Њ "%Ђ"

0џџџџџџџџџ
 
,__inference_generator_l2_layer_call_fn_71845O#$/Ђ,
%Ђ"
 
inputsџџџџџџџџџ

Њ "џџџџџџџџџЇ
G__inference_generator_l3_layer_call_and_return_conditional_losses_71856\%&/Ђ,
%Ђ"
 
inputsџџџџџџџџџ
Њ "%Ђ"

0џџџџџџџџџ
 
,__inference_generator_l3_layer_call_fn_71865O%&/Ђ,
%Ђ"
 
inputsџџџџџџџџџ
Њ "џџџџџџџџџЛ
D__inference_generator_layer_call_and_return_conditional_losses_70837s!"#$%&'(@Ђ=
6Ђ3
)&
generator_inputџџџџџџџџџ
p

 
Њ "%Ђ"

0џџџџџџџџџ
 Л
D__inference_generator_layer_call_and_return_conditional_losses_70861s!"#$%&'(@Ђ=
6Ђ3
)&
generator_inputџџџџџџџџџ
p 

 
Њ "%Ђ"

0џџџџџџџџџ
 В
D__inference_generator_layer_call_and_return_conditional_losses_71625j!"#$%&'(7Ђ4
-Ђ*
 
inputsџџџџџџџџџ
p

 
Њ "%Ђ"

0џџџџџџџџџ
 В
D__inference_generator_layer_call_and_return_conditional_losses_71657j!"#$%&'(7Ђ4
-Ђ*
 
inputsџџџџџџџџџ
p 

 
Њ "%Ђ"

0џџџџџџџџџ
 
)__inference_generator_layer_call_fn_70907f!"#$%&'(@Ђ=
6Ђ3
)&
generator_inputџџџџџџџџџ
p

 
Њ "џџџџџџџџџ
)__inference_generator_layer_call_fn_70952f!"#$%&'(@Ђ=
6Ђ3
)&
generator_inputџџџџџџџџџ
p 

 
Њ "џџџџџџџџџ
)__inference_generator_layer_call_fn_71678]!"#$%&'(7Ђ4
-Ђ*
 
inputsџџџџџџџџџ
p

 
Њ "џџџџџџџџџ
)__inference_generator_layer_call_fn_71699]!"#$%&'(7Ђ4
-Ђ*
 
inputsџџџџџџџџџ
p 

 
Њ "џџџџџџџџџЋ
K__inference_generator_output_layer_call_and_return_conditional_losses_71876\'(/Ђ,
%Ђ"
 
inputsџџџџџџџџџ
Њ "%Ђ"

0џџџџџџџџџ
 
0__inference_generator_output_layer_call_fn_71885O'(/Ђ,
%Ђ"
 
inputsџџџџџџџџџ
Њ "џџџџџџџџџЃ
#__inference_signature_wrapper_71381|)*+,-./0;Ђ8
Ђ 
1Њ.
,
input_1!
input_1џџџџџџџџџ"3Њ0
.
output_1"
output_1џџџџџџџџџ