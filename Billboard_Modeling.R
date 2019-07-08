############################################################################
## Data :  21/03/2018)                                                    ##
## Autor : Pedro Henrique Sodré Puntel                                    ##
## Email : pedro.puntel@gmail.com                                         ##
## Instituição : Escola Nacional de Ciencias Estatísticas - ENCE IBGE     ##
## Disciplina : Projeto de Iniciaçãoo Científica - PIBIC CNPq 2018/2019   ##
## Professor : Gustavo Ferreira                                           ##
## Tema : Aplicação dos Modelos para redes sociais aos dados da Billboard ##
############################################################################

# Encoding default para este script é UTF-8 

# Pacotes utilizados
library(dplyr)
library(plotly)
library(network)
library(MCMCpack)
library(latentnet)
library(compiler)

# Perfomance Setup
cl = makeCluster(detectCores(), type = "PSOCK")
registerDoSEQ()

# Importa a sociomatriz
path = "C:\\Users\\pedro\\Desktop\\R\\PIBIC 2018\\Analyses\\BillBoard\\Database\\15062019_BillboardTop20_SocialMatrix.csv"
Social_Matrix = rio::import(path) %>% as.matrix()
rownames(Social_Matrix) = colnames(Social_Matrix)
rm(path)

# Redução da sociomatriz - "100 artistas que mais ocuparam o Top 20 da Billboard"
Social_Matrix = Social_Matrix[order(rowSums(Social_Matrix), decreasing = T), order(colSums(Social_Matrix), decreasing = T)]
Social_Matrix = Social_Matrix[1:100,1:100]

#############################
## Estatísticas Descritivas : 
############################# 

# Heatmap da sociomatriz
plot_ly(z = Social_Matrix, x = colnames(Social_Matrix), y = colnames(Social_Matrix),
        type = "heatmap", zauto = F, zmin = 0, zmax = 64, colors = c("Yellow","Orange","Red")) %>%
    layout(title = "Heatmap Sociomatriz",
           xaxis = list(categoryorder = "array", categoryarray = colnames(Social_Matrix)),
           yaxis = list(categoryorder = "array", categoryarray = colnames(Social_Matrix))) %>%
    colorbar(title = "Meses no Top 20")
  
# Artistas que mais ocuparam o Top 20 da Billboard
aux.vec <- diag(Social_Matrix) %>% as.numeric()
aux.vec <- aux.vec[order(aux.vec, decreasing = T)]
aux.df <- cbind(colnames(Social_Matrix), aux.vec) %>% data.table::as.data.table()
plot_ly(aux.df, x = colnames(Social_Matrix), y = aux.vec, type = 'bar',
        marker = list(color = "Purple", line = list(color = "Purple", width = 2.5))) %>%
    layout(title = "Número de meses ocupados no Top 20 da Billboard - (65 no total)",
           xaxis = list(categoryorder = "array", categoryarray = colnames(Social_Matrix)),
            yaxis = list(title = "Meses no Top 20", categoryorder = "array", categoryarray = aux.vec))

# Limpando memória
rm(aux.df, aux.vec)

############## Análises ###############
# . Heatmap mostra claramente que poucos são os artistas que ocuparam o Top 20 da Billboard por muito tempo
# . Destaque para os Top 5 artistas : Drake, Ed Sheeran, Taylor Swift, Ariana Grande e The Weekend
# . Gráfico de barras reforça o observado no Heatmap.
  
##########################################
## I - Modelos para abertura de Conexões :
########################################## 

# [03/06/19] :
# Chegou-se a conclusão de que este modelo deve ser desconsiderado uma vez que a sua
# interpretação poder ser facilmente substituida pelos gráficos da sessão de Estatística Descrtiva.

# Parâmetros estimados pelo modelo podem ser interpretados como o grau de extroversão/introversão
# de cada artista, ou seja, o quão propício um determinado artista é de compartilhar o Top 20 da
# da Billboard com um outro artista qualquer.
#
# A estimação dos parâmetros é feita por meio de Modelos Lineares Generalizados (MLG's). No R, antes
# mesmo de chamar a função glm()  aos dados, é necessári a criação de uma "Matriz de Delineameto", cujo
# propóstio é o de servir como uma "plataforma" para os parâmetros. Por exemplo, no nosso caso com a 
# sociomatriz Billboard, existem Comb(ncol, 2) parâmetros a serem estimados. Assim, a Matriz de Delineamento
# basicamente informa as funçõess lm() e glm() do R quantos parâmetros devem ser estimados.
#
#  . A dimensão da Matriz de Delineamento é sempre Comb(ncol,2) x ncol
#  . Suas linhas correspondem sempre ao triângulo superior da sociomatriz em questão
#
# Existe ainda a Matriz de Delineamento que leva em consideração o intercepto. Porém, nosso modelo não
# considera um intercepto e então esta versão da matriz não foi considerada.

# Função que implementa o Modelo de Abertura para Conexões
Compute.Modelo_Abertura_Conexoes = function(socio.mat, k.artists, n.months) {
  
  # Sub-rotina responsável pela construção da Matriz de Delineamento
  Build.Design_Matrix = function(k.artists) {
    
    n = k.artists*(k.artists-1)/2
    X =  matrix(0,n,k.artists)
    aux = matrix(0,n,2)
    cont = 1
    
    for (i in 1:(k.artists-1)){
      for (j in (i+1):k.artists){
        aux[cont,] = c(i,j)
        cont = cont+1
      }
    }
    
    aux
    
    for (i in 1:n){
      X[i,aux[i,1]] = 1
      X[i,aux[i,2]] = 1
    }
    
    return(X)
  }
  
  # Sub-rotina responsável pela construção Matriz de Mínimos
  Build.Minimum_Matrix = function(socio.mat) {
    
    # IDefine a estrutura da Matriz de Mínimos associada
    min.mat = matrix(NA, ncol = ncol(socio.mat), nrow = nrow(socio.mat))
    
    # Diagonal da Sociomatriz
    m.diag = diag(socio.mat)
    
    ## Constrói a Matriz de Mínimos
    for (i in 1:ncol(socio.mat)) {
      for (j in 1:nrow(socio.mat)) {
        min.mat[i,j] = min(m.diag[i], m.diag[j])
      }
    }
    
    # Nomeia as linhas e as Colunas da Matriz de Mínimos
    colnames(min.mat) = colnames(socio.mat)
    rownames(min.mat) = rownames(socio.mat)
    
    # Retorna a Matriz de Mínimos
    return(min.mat)
    
  }
  
  # Sub-rotina que implementa o modelo assumindo : Y_ij ~ Po(X_ij)
  Poisson.Model = function(socio.mat) {
    
    # Vetoriza as entradas o triângulo inferior da Sociomatriz
    m.vec = socio.mat[lower.tri(socio.mat, diag = F)]
    
    # Constrói a matriz de delineamento 
    design.mat = Build.Design_Matrix(ncol(socio.mat))
    
    # Estimação
    Glm.Poisson = glm(m.vec ~ design.mat - 1, family = "poisson")
    
    # Tabela auxiliar
    alphas = Glm.Poisson$coefficients
    aux.mat = cbind(colnames(socio.mat), alphas) %>% data.table::as.data.table()
    ## aux.mat <- aux.mat[order(-alphas),]
    rownames(aux.mat) = as.character(seq(1,nrow(socio.mat),1))
    colnames(aux.mat) = c("Names","Alphas")
    
    # Gráfico dos parâmetros estimados pelo modelo
    Glm.Poisson_Plot = plot_ly(aux.mat, x = colnames(socio.mat), y = alphas, type = 'bar',
            marker = list(color = "Green", line = list(color = "Green", width = 2.5))) %>%
      layout(title = "Modelo de Abertura de Conexões - MLG Poisson",
             xaxis = list(title = "Artistas", categoryorder = "array", categoryarray = aux.mat$Names),
             yaxis = list(title = "Grau de Abertura", categoryorder = "array", categoryarray = aux.mat$Alphas))
    
    # Retorna as informaçãoes
    list("Model" = Glm.Poisson, "Plot" = Glm.Poisson_Plot) %>% return()
    
  }
  
  # Sub-rotina que implementa o modelo assumindo : Y_ij ~ Bin(N_ij, p_ij)
  Binomial.Min_Model = function(socio.mat) {
    
    ## Modificação do modelo com a introdução da matriz de mínimos
    ##
    ##          Y_ij | P_ij ~ Binomial (min(n_i,n_j), p_ij)
    ##
    ## Onde min(n_i,n_j) pode ser interpretado como uma ponderação
    ## acerca do númeorde vezes que o i-ésimo e o j-ésimo artista
    ## da rede ocuparam o Top 20. A motivação por trás do mínimo é
    ## simples : "Se i ocupou o Top 20 muito mais vezes do que j,
    ## entãoo a relação entre i e j deve ser pautada somente pelo
    ## número de vezes em que j ocupou o Top 20".
    ##
    ## Do contrário, obviamente, i terá sempre um grau de abertura maior.
    ## Este era justamente o vício introduzido pelo modelo poisson.
    ##
    ## Note que a implementação do Modelo Linear Generalizado no R para a 
    ## família Binomial é diferente daquela para a família Poisson : A função
    ## glm() recebe, neste caso, uma matriz de dimensão C(K,2) x 2 cujas colunas
    ## representam, 'núero de sucessos' e  o 'número de fracassos', respectivamente.
    
    ## Vetor de sucessos
    p.vec <- socio.mat[lower.tri(socio.mat)]
    
    ## Constróia matriz de mínimos
    min.mat <- Build.Minimum_Matrix(socio.mat)
    
    ## Vetor de fracassos
    q.vec <- min.mat[lower.tri(min.mat)] - p.vec
    
    ## Constrói a Matriz de Delineamento
    design.mat <- Build.Design_Matrix(k.artists)
    
    ## Estimação
    Glm.Binomial <- glm(cbind(p.vec,q.vec) ~ design.mat - 1, family = "binomial")

    ## Tabela auxiliar
    Alphas <- Glm.Binomial$coefficients
    mat.aux <- cbind(colnames(socio.mat), Alphas) %>% data.table::as.data.table()
    mat.aux <- mat.aux[order(-Alphas),]
    rownames(mat.aux) <- as.character(seq(1,nrow(socio.mat),1))
    colnames(mat.aux) <- c("Names","Alphas")
    
    ## Gráfico dos parâmetros estimados pelo modelo
    Glm.Binomial_Plot <- plot_ly(mat.aux, x = colnames(socio.mat), y = Alphas, type = 'bar',
            marker = list(color = 'Light Blue',line = list(color = 'Light Blue', width = 2.5))) %>%
      layout(title = "Modelo de Abertura de Conexões - MLG Binomial 1",
             xaxis = list(title = "Artistas", categoryorder = "array", categoryarray = mat.aux$Names),
             yaxis = list(title = "Grau de Abertura", categoryorder = "array", categoryarray = mat.aux$Alphas))
    
    ## Retorna as informações
    list("Model" = Glm.Binomial, "Plot" = Glm.Binomial_Plot) %>% return()
    
  }
  
  ## Sub-rotina que implementa o modelo assumindo : Y_ij ~ Bin(n.months, p_ij)
  Binomial.Tot_Model <- function(socio.mat, n.months) {
    
    ## Dado que o número máximo que um determinado artista poderia assumir 
    ## na nossa rede é justamente o número meses scrapeados, o que aconteceria
    ## com as estimativas do modelos se considerássemos estes valores como sendo
    ## o nosso novo vetor de sucessos no modelo Binomial ?
    ##
    ## Esta suposição parte do princípio que todo artista detêm o do mesmo potencial
    ## para ocupar o Top 20 da Billboard. Algo que minimamente contestável.
    
    ## Vetor de sucessos
    p.vec <- socio.mat[lower.tri(socio.mat)]
    
    ## Vetor de fracassos
    q.vec <- rep(n.months,length(p.vec)) - p.vec
    
    ## Constrói  Matriz de Delineamento
    design.mat <- Build.Design_Matrix(k.artists)
    
    ## Estimação dos Parâmetros do Modelo Linear Generalizado
    Glm.Binomial_2 <- glm(cbind(p.vec, q.vec) ~ design.mat - 1, family = "binomial")
    
    ## Data Table auxiliar
    alphas <- Glm.Binomial_2$coefficients
    mat.aux <- cbind(colnames(socio.mat), alphas) %>% data.table::as.data.table()
    mat.aux <- mat.aux[order(-alphas),]
    rownames(mat.aux) <- as.character(seq(1,nrow(socio.mat),1))
    colnames(mat.aux) <- c("Names","Alphas")
    
    ## Gráfico dos parâmetros estimados pelo modelo
    Glm.Binomial_Plot_2 <- plot_ly(mat.aux, x = colnames(socio.mat), y = alphas, type = 'bar',
            marker = list(color = 'Purple',line = list(color = 'Purple', width = 2.5))) %>%
      layout(title = "Modelo de Abertura de Conexões - MLG Binomial 2",
             xaxis = list(title = "Artistas", categoryorder = "array", categoryarray = mat.aux$Names),
             yaxis = list(title = "Grau de Abertura", categoryorder = "array", categoryarray = mat.aux$Alphas))
    
    ## Retorna as informações
    list("Model" = Glm.Binomial_2, "Plot" = Glm.Binomial_Plot_2)
    
  }
  
  ## Retorna os resultados
  list("Poisson" = Poisson.Model(socio.mat),
       "Min.Binom" = Binomial.Min_Model(socio.mat),
       "Tot.Binom" = Binomial.Tot_Model(socio.mat, n.months)) %>% return()
  
}
Compute.Modelo_Abertura_Conexoes <- cmpfun(Compute.Modelo_Abertura_Conexoes)
MLG <- Compute.Modelo_Abertura_Conexoes(Social_Matrix, k.artists = 100, n.months = 64)
which(as.numeric(MLG$Poisson$Model$coefficients) > 0) %>% length()
MLG$Poisson$Model$residuals %>% as.numeric() %>% plot()
MLG$Poisson$Model$converged
MLG$Poisson$Plot

##############  Resultados ###############
## . Modelo Poisson aparenta ser o mais relevante dentre os três
## . Modelo Poisson foi o único que convergiu dos três
## . Modelo Poisson classifica 46 artistas como extrovertidos e 54 como introvertidos
## . Não se observou normalidade em nenhum dos resíduos dos três modelos
## . Uma simples ordenação dos artistas pela soma das suas linhas daria a mesma conclusão

#################################################
## II - Escalonamento Multidimensional Clássico :
################################################# 

## O Escalonamento Multidimensional Clássico supõem existência de uma matriz de distâncias
## D[n x n], cujas entradas d_ij quantificam a distância entre o elemento i e o elemento j
## sob o contexto dos dados originais.
##
## O objetivo do método é então encontrar uma configuração de pontos no espaço p-dimensional,
## de tal forma que as coordenadas dos pontos, ao longo da dimensão p, produzam uma matriz de
## distâcias Euclidianas, cujos elementos estejam mais próximos possível dos de D.
##
## A justificativa para a implementação deste método é imediata :
## 
##       Gostaríamos de verificar a similiarede/dissimilarides dos artistas da rede
##          quanto ao seu padrão de conexões ou seja, quais seriam os artistas que
##                        mais assemelham-se em termos de perfil ? 
##
## Tão imediata quanto a sua justificativa é a sua interpretação :
##
##   Um grupo de dois ou mais artistas podem ser considerados semelhantes se a distância 
##      entre estes no espaço euclidano p-dimensianal (p <= n) pequena e dissimilares
##                                  caso contrário.

## [03/06/19] : Chegou-se a conclusão de que o Escalonamento Multidimensional, quando 
## aplicado para com uma matriz de adjacência valorada, acaba introduzindo o seguinte
## viés : atores de uma rede cujo padrão de conexões é similar, acabam ficando muito
## próximos uns dos outros por mais que estes mesmos, entre si, não estabelecam nenhum
## tipo de relação. No nosso contexto, por exemplo, pode-se verificar este efeito com
## os artistas "não-populares" que ficam todos próximos uns dos outros (aglomerado),
## mesmo que entre eles, não existam conexões. A saída encontrada foi realizar o 
## Escalonamento Multidimensional com uma "Matriz de Correlações" cujas entradas
## quantificam o Coeficiente de Correlação de Pearson entre os artistas i e j.

## Função que implementa o Escalonamento Multidimensional
Compute.Multidimensional_Scaling <- function(socio.mat, use.risk_mat = F) {
  
  ## Sub-rotina responsável pela construção Matriz de Mínimos
  Build.Minimum_Matrix <- function(socio.mat) {
    
    ## Define a estrutura da Matriz de Mínimos associada
    min.mat <- matrix(NA, ncol = ncol(socio.mat), nrow = nrow(socio.mat))
    
    ## Diagonal da Sociomatriz
    m.diag <- diag(socio.mat)
    
    ## Constrói a Matriz de Mínimos
    for (i in 1:ncol(socio.mat)) {
      for (j in 1:nrow(socio.mat)) {
        min.mat[i,j] <- min(m.diag[i], m.diag[j])
      }
    }
    
    ## Nomeia as linhas e as Colunas da Matriz de Mínimos
    colnames(min.mat) <- colnames(socio.mat)
    rownames(min.mat) <- rownames(socio.mat)
    
    ## Retorna a Matriz de Mínimos
    return(min.mat)
    
  }
  
  ## Implementa o Escalonamento Multidimensional Clássico
  if (isFALSE(use.risk_mat)) {
    
    ## Computa o Escalonamento Multidimensional com a matriz de correlações
    Mds <- dist(cor(socio.mat), method = "euclidean", diag = T, upper = T) %>% cmdscale(k = 2) %>% as.matrix()
    
    ## Nomeia as colunas do Data Table resultate para melhor manipulação
    colnames(Mds) <- c("x.coord","y.coord")
    
    ## Visualização da rede em forma de Grafo
    ## [!] Divide-se a Socio-Matriz por 10 para reduzir o número de arestas plotadas
    edges <- igraph::graph.adjacency(as.matrix(round(socio.mat/10,0)), mode = "undirected", diag = F)
    Mds.plot <- igraph::tkplot(edges, vertex.size = 30, edge.width = 1, vertex.label = colnames(socio.mat),
     layout = Mds, vertex.color = I("Dark Orange"), vertex.label.cex = 1, vertex.label.color = "Black")

  }

  ## Implementa o Escalonamento Multidimensional com a Matriz de Riscos
  else {
    
    ## Suponha que uma epidemia atingiu 4 cidades distintas e gostaríamos de quantificar
    ## qual destas foi a mais tingidda/contaminda. É atural que de prontidão, a nossa
    ## solução seja contar o número de contágios da doença que ocorreram em cada uma
    ## das cidades e assim, verificar qual delas foi a mais impactada.
    ##
    ## Porém, existe um problema com essa abordagem: estamos desconsiderando a população
    ## de cada cidade para efeito de comparação. Isto pode ser potencialmente problemático
    ## uma vez que 50 casos de contágio em um população de 100 pessoas é muito mais preocupante
    ## do que 500 casos em uma população de 100000 pessoas. 
    ##
    ## É justamente neste sentido que convêm realizarmos o Escalonamento Multidimensional não
    ## diretamente ao nosso conjunto de dados mas sim com uma "Matriz de Riscos". Seu procediemnto
    ## de montagem é descrito a seguir :
    ##
    ##  Seja A(nxn) a sociomatriz de interesse. De forma a entender proporcionalmente o quão
    ##  forte ou fraca é a conexão entre dados dois atores da rede, montemos uma nova matriz
    ##  P(nxn) onde cada entrada p_ij quantifica o quão acima/abaixo da média (ou do esperado)
    ##  se dá a relação entre tais atores i e j. 
    ##    
    ##  Para tal, definemos uma taxa de normalização alpha, dada pela razão entre a soma das
    ##  entradas da sociomatriz A, pela a soma das quantidades as os tais atores poderiam,
    ##  supostamente, estabelecer uns com os outros.
    ##    
    ##  Assim, preenchemos por fim a matriz P onde cada entrada p_ij seja então a razão entre
    ##  as entradas da sociomatriz (valor real) pelos valores esperados para tais conexões
    ##  (valor teórico).
    ##
    ## Interpretação : "Quão acima da média se dá a relação entre os atores i e j da rede ?"
    
    ## Computa a soma das entradas do triângulo inferior da sociomatriz
    sum.conex <- socio.mat[lower.tri(socio.mat, diag = F)] %>% sum()
    
    ## Constrói a Matriz de  Mínimos
    min.mat <- Build.Minimum_Matrix(socio.mat)
    
    ## Computa a soma das entradas do triângulo inferior da matriz de mínimos
    sum.min_conex <- min.mat[lower.tri(min.mat, diag = F)] %>% sum()
    
    ## Calculo da Taxa de Normalização
    ## . Valor entre 0 e 1
    ## . --- > 1 indicam que estabeleceram-se mais conexões do que o esperado
    ## . <--- 1 indicam que estabeleceram-se menos conexões do que o esperado
    taxa.norm <- sum.conex / sum.min_conex
    
    ## Montagem da Matriz de Riscos
    risk.mat <- taxa.norm * min.mat
    colnames(risk.mat) <- colnames(socio.mat)
    rownames(risk.mat) <- rownames(socio.mat)
    
    ## Computa o Escalonamento Multidimensional
    Mds <- dist(risk.mat, method = "euclidean", diag = F, upper = T) %>%
      cmdscale(k = 2) %>% as.matrix()
    
    ## Nomeia as colunas do Data Table resultate para melhor manipulação
    colnames(Mds) <- c("x.coord","y.coord")
    
    ## Visualização em forma de Grafo
    edges <- igraph::graph.adjacency(as.matrix(round(socio.mat/10,0)), mode = "undirected", diag = F)
    Mds.plot <- igraph::tkplot(edges, vertex.size = 30, edge.width = 1, vertex.label = colnames(socio.mat),
     layout = Mds, vertex.color = I("Dark Orange"), vertex.label.cex = 1, vertex.label.color = "Black")
    
  }
  
  ## Retorna os resultados
  list("Mds.par" = Mds) %>% return()
  
}
Compute.Multidimensional_Scaling <- cmpfun(Compute.Multidimensional_Scaling)
MDS <- Compute.Multidimensional_Scaling(Social_Matrix[1:100,1:100], use.risk_mat = F)

############## Resultados ################
## . Padrão de conexões ainda é conservado após a redução da Socio-Matriz.
## . Descartou-se o Escalonamento Multdimensional feito com a matriz de riscos.
## . O Escalonamento Multidimensional feito com a Socio-Matriz de correlações
## realmente contorna o vício mencionado, expondo relações interessantes entre
## os artistas.Nesse sentido, convêm destacar a relação do artista Drake com
## o artista Pentatonix, por exemplo. Note que apesar do primeiro ser o mais
## popular, ou seja, apesar do Drake ter alta correlação com todos os demais
## artistas, artistas como Pentatonix, Metallica e Rae Sremmurd são aqueles
## mais próximos do Drake pelo simples de fato de que nas poucas vezes em
## que estes estiveram no Top 20 da Billboard, o Drake certamente ocupava
## o Top 20 também.

########################################
## III - Modelo de Distâncias Latentes :
######################################## 

## Os modelos de Espaços Latentes propõem uma classe de modelos nos quais a probabilidade
## de uma relação entre os atores de uma rede é uma função não-linear das suas respectivas
## posições em um espaço intitulado "espaço social". Dentre as possíves interpretações para
## o referido espaço, existe aquela de um espaço não-observacional que retrata o padrão de
## conectividade dos atores da rede como base nas suas posições. Intuitivamente, atores mais
## próximos são mais próprios de estabelecerem uma relação do que aqueles mais afastados.
## 
## Estatisticamente, o modelo caracteriza-se por assumir que cada vértice possui uma posição
## desconhecida no espaço social, e que uma relação entre dois vértices são condicionalmente
## independentes dadas estas suas respectivas posições. Assim, a probabilidade de uma conexão
## entre estes indivíduos é modelada no seguinte formato: 
##
##                                  logit(pij) ~ Q - |Zi-Zj|
##
##               * Q = nível de atividade social de cada ator (covariável) 
##               * |Zi-Zj| = distância entre os vertices i e j no espaço social
##
## 'Q' pode ser pensado também como o grau de extroversão do i-ésimo ator da rede, de forma
## que a distância entre um ator aos demais seja "ponderada" também pelo quão propício tal ator
## é de estabelecer uma conexão.O
##
## Cada 'Zi' e 'Zj' são vetores de dimensão 2 (coordenadas no plano bi-dimensional).
##
## As distâncias |Zi-Zj| entre os atores são invariantes quanto a reflexão, rotação ou deslo-
## camento dos eixos coordenados. Como consequência, é comum fixar 3 dos parâmetros para fins de
## estimação.
##
## Valores iniciais candidatos os posições dos atores no espaço latente são aqueles proveni-
## entes do escalonamento multidimensional. Tais valores serão posteriormente otimizados pela
## função maximizadora não-linear do R nlminb().

## [03/06/19] : Por uma mera questão de consistência com a mudança introduzida no Escalonamento
## Multidimensional, optou-se por realizar a estimação dos parâmetros utilizando a Socio-Matriz
## de Correlação (tanto nos métodos de Máxima Verossimilhança como nos métodos Bayesianos).

## Função que implementa o Modelo de Distâncias Latentes 
Compute.Latent_Distances_Model <- function(socio.mat, Is.direct, Is.Weighted, Has.Loops, Bayesian.Est) {

  ## Sub-rotina que implementa a estimação por Máxima Verossimilhança
  MLH.Estimation <- function(socio.mat) {
    
    ##                    --> [Yij | Lij] ~ Poisson(Lij)
    ##
    ##                    --> Lij = Exp(Q - |Zi - Zj|)
    ##
    ##  --> L(Q, Zi, Zj | Yij) = Sum(i<j) { Yij*ln(Lij) - Lij - ln(Yij!) }
    ##
    ## . Lij = Exp(Q - |Zi-Zj|) é sempre o valor do parâmetro Lambda entre os atores i e j da rede
    ## . Sum(i<j) percorre o triângulo superior da sociomatriz Y
    ## . O parâmetro Q (theta) é como se fosse uma distância média entre todos os atores da rede.
    ## . Para |Zi-Zj| grande, a probabilidade de uma relação entre os atores i e j é pequena.
    ## . A última entrada do vetor de parâmetros é sempre fixa e corresponde ao parâmetro 'Q'
    ## . sqrt(...) remete a fórmula da distância euclidianda entre dois pontos
    ## . Fixa-se sempre os valores de Z11,Z21,Z31 (invariância quanto a rotação, refelxão...)
    
    ## Constrói o vetor de parâmetros a partir do Escalonamento Multidimensional com a Sociomatriz de Correlações
    par.vec <- dist(cor(socio.mat), method = "euclidean", diag = F, upper = T) %>%  cmdscale(k = 2) %>% as.vector()
    par.vec[(length(par.vec)+1)] <- 0
    
    ## Função de Log-Verossimilhança
    Poisson.LLH <- function(x) { 
      
      ## Fixa os três primeiros parâmetros (exceto o Q)
      x[1] <- par.vec[1]
      x[2] <- par.vec[2]
      x[3] <- par.vec[3]
      
      ## Número de atores na rede 
      K <- ncol(socio.mat)
      
      ## Variável auxiliar que guarda o valor da função calculada para cada par de atores da rede
      fun.val <- NULL
      
      ## Itera por todo o triângulo superior da Sociomatriz (exceto diagonal)
      for(i in 1:(K-1)) { 
        for(j in (i+1):K) { 
          
          ## Cálculo do parâmetro Lambda associado a cada par de atores da rede
          lambda <- exp(x[2*K+1]-sqrt((x[i]-x[j])^2+(x[i+K]-x[j+K])^2))
          
          ## Guarda o valor da função em cada iteração
          fun.val <- c(fun.val, socio.mat[i,j]*log(lambda) - lambda - log(factorial(socio.mat[i,j])))
          
        }
      }
      
      ## Retorna a soma dos valores da função de Log-Verossimilhança
      return((-1)*sum(fun.val))
      
    }
    
    ## Maximização da função de Log-Verossimilhança
    Poisson.MLE <- nlminb(start = par.vec, objective = Poisson.LLH)
    
    ## Tabela com os parâmetros estimados
    MLE.coord <- cbind(Poisson.MLE$par[1:ncol(socio.mat)], Poisson.MLE$par[(ncol(socio.mat)+1):(2*ncol(socio.mat))])
    rownames(MLE.coord) <- colnames(socio.mat)
    
    ## Visuzalização dos parâmetros estimados em forma de Grafo
    edges <- igraph::graph.adjacency(as.matrix(round(socio.mat/10,0)), mode = "undirected", diag = F, weighted = T)
    igraph::tkplot(edges, vertex.size = 30, edge.width = 1, vertex.label = colnames(socio.mat),
            layout = MLE.coord, vertex.color = I("Light Blue"), vertex.label.cex = 1, vertex.label.color = "Black")
    
    ## Retorna as informaçõees
    list("Has.Converged" = Poisson.MLE$message, "Fun.Maximum" = Poisson.MLE$objective,
         "Fun.Eval" = Poisson.MLE$evaluations, "MLE.Par" = MLE.coord) %>% return()
    
  }
  
  ## Sub-rotina que implementa a estimação por métodos Bayesianos
  Bayesian.Estimation <- function(socio.mat, Is.direct, Is.Weighted, Has.Loops) {
    
    ## Método de simulação estocástica cujo objetivo é obter amostras da distribuição
    ## "a posteriori" quando as mesmas não possuem um forma fechada, ou seja, quando não
    ## pertencem à nenhuma distribuição conhecida. De forma simples, o método busca estimar
    ## a função de distribuição "a posteriori" possibilitando assim inferência acerca dos
    ## parâmetros de interesse.
    ##
    ## O objetivo do método de MCMC é construir uma cadeia de Markov cuja distribuição limite
    ## seja igual a distribuição de interesse. O algortimo realiza um número finito de 
    ## sucessivas de amostagens desta cadeia, de forma a esperar-se que a sua distribuição
    ## convirja para a tal distribuição "a posteriori". 
    ##
    ## Dentre os métodos mais utlizados para a construção desta cadeia de Markov, destacam-se
    ## o "Amostrador de Gibbs" e o método de "Metropolis-Hastings" :
    ##
    ## O método Amostrador de Gibbs é um método iterativo de amostragem de uma cadeia de Markov,
    ## cujas probabilidades de transição são definidas de acordo com a distribuição condicional
    ## completa dos parâmetros | conjunto de dados observados.
    ##
    ## O método Metropolis-Hastings é um método que extrai amostras de qualquer distribuição de
    ## probabilidade P(x), dado o valor de uma função f(x) que seja proporcional é densidade P.
    ## O algortimo funciona de forma iterativa, amostrando da função proporcional f(x), de forma
    ## que o próximo valor amostrado seja exclusivamente dependente do seu anterior, formando 
    ## assim uma sequência de amostras que concomitam em uma cadeia de Markov. A diferença,
    ## porém, reside em uma espécie de "filtro" que avalia a plausabilidade dos valores
    ## amostrados, de forma que onde em cada iteração, o algortimo aceita/rejeita determinada
    ## amostra com certa probabilidade. Caso aceita, seu valor será utilizado na próxima
    ## iteração e caso contrário, descarta-se este mantendo aquele da iteração anterior.
    ##
    ## Assim, os processos de estimação descritos acima serão aplicados no R com o auxílio dos
    ## pacotes 'network', 'latentnet' e 'MCMCpack' como segue :
    
    ## Objeto tipo 'network' a ser manipulado pela função ergmm()
    net.aux <- network(socio.mat, directed = Is.direct, loops = Has.Loops, matrix.type = "adjacency")
    
    ## Atribui pesos as arestas da rede (caso esta seja de fato valorada)
    ## . http://doogan.us/netdata.html
    if (isTRUE(Is.Weighted)) {
      set.edge.value(net.aux, attrname = "weight", value = socio.mat)
    }

    ## Exponential Random Graph Mixed Models (ERGMM's)
    ## 
    ## . pmode = Posterior mode (MAP)
    ## . mcmc = Monte Carlo Markov Chain
    ## . mkl = Kullback-Leiber weighted likelihood
    ## . mle = Maximum Likelihood Estimation
    ## . procrustes = Procrustean Analysis
    ##
    ## Dentre os ERGMM's aplicados, nos interessa somente as estimações referente é moda
    ## a posteriori, Máxima Verossimilhança e Procrustes. Uma breve explicação destes
    ## modelos é feita a seguir (exceto sobre a Máxima Verossimilhança que já sabemos) :
    ##
    ## Uma estimativa de Máxima Densidade a Posteriori (MAP) é uma estimativa de uma
    ## quantidade desconhecida, que é igual a moda da distribuição posterior. O MAP pode ser
    ## usado para obter uma estimativa pontual de uma quantidade não observada com base em
    ## dados empíricos. Esta é relacionada ao método de estimação por máxima
    ## verossimilhança, mas emprega um objetivo de otimização que incorpora uma distribuição
    ## a priori (esta fornece informações adicionais através de um evento prévio) sobre a
    ## quantidade que se deseja estimar.

    ## Procrustes Analysis
    ##
    ## Distances between a set of points in Euclidean space are invariant under rotation,
    ## reflection, and translation. Therefore, for each matrix of latent positions Z, there
    ## is an infinite number of other positions giving the same log-likelihood.
    ##
    ## A confidence region that includes two equivalent positions Z1 and Z2 is, in a sense,
    ## overestimating the variability in the unknown latente positions, because these are
    ## identical for Z1 and Z2.
    ##
    ## Fortunately, this problem can be resolved by perfoming inference on a seta that we call 
    ## "equivalence classes of latent positions" : Let [Z] be the set of positions equivalent
    ## to matrix of latent positions Z under the previous mentioned operations. For each [Z],
    ## there is one set of distances between the nodes. We oftne refet to this class of positions
    ## as a "configuration model".
    ## 
    ## Therefore, We make inference on these configurations by performing inference on particular
    ## elements that are comparable across configurations. For example, given an configuration set
    ## [Z], we select for inference the configuration Z* that minimizes the sum of squared positional
    ## difference. It turns out that Z* ends up beign a "Procrustean" transformation of Z, therefore
    ## being the element of the set [Z] which is closest to the the original "fixed" set.
    ##
    ## We typically take for this orginal "fixed" set the MLE estimates of the latent positions.
    
    ## [03/06/09] : Valores iniciais provenientes do Escalonamento Multidimensional com a matriz de correlações
    par.vec <- dist(cor(socio.mat), method = "euclidean", diag = F, upper = T) %>%  cmdscale(k = 2) %>% as.vector()
    par.vec[(length(par.vec)+1)] <- 0
    
    ## Fits the Model
    Bayesian.Fit <- ergmm(net.aux ~ euclidean(d = 2), family = "Poisson", user.start = list(par.vec),
                          tofit = c("pmode","mcmc","mkl","mkl.mbc","mle","procrustes","klswitch"))
    
    ## Retorna as informações
    list("Results" = Bayesian.Fit) %>% return()
    
  }
  
  ## Retorna as informações
  if (isFALSE(Bayesian.Est)) {
    list("Results" = MLH.Estimation(socio.mat)) %>% return()
  } else {
    list("Results" = Bayesian.Estimation(socio.mat, Is.direct, Is.Weighted, Has.Loops)) %>% return()
  }
  
}
Compute.Latent_Distances_Model <- cmpfun(Compute.Latent_Distances_Model)

############## Resultados (Máxima Verossimilhança) #############################
## . Infelizmente, só foi possível maximizar a função de Log-Verossimilhança com 25 artistas
## . Mesmo aumentando o número de iterações para a função nlminb(), não houve diferença
## . Disposição dos artistas distinta daquela sugerida pela Escalonamento Multidimensional
## . É difícil identificar um padrão de agrupamento com apenas 25 artistas, porém, mesmo assim 
## é possível identificar os artistas mais populares aglomerados no centro da rede enquanto os
## demais (menos populares) vão se posicionando ao reder destes.
Latent.MLE <- Compute.Latent_Distances_Model(Social_Matrix[1:25,1:25],F,T,T,F)
Latent.MLE$Results$Has.Converged
Latent.MLE$Results$Fun.Maximum
Latent.MLE$Results$MLE.Par

############## Resultados (Bayesianos) #############################
## . Métodos Bayesianos já são capazes de lidar com a Sociomatriz por completo
## . Valores passados como parâmetros iniciais pouco interferem na disposição final dos artistas
## . par.vec = MDS ou par.vec = cor.MDS resultam em grafos (es estimações) parecidissímas
## . Parâmetros estimados por MV diferem numericamente mas não graficamente daqueles por MDP
## . Disposição dos artistas totalmente diferente daquela sugerida pelo Escalonamento Multidimensional
Latent.Bayesian <- Compute.Latent_Distances_Model(Social_Matrix[1:25,1:25],F,T,F,T)
Latent.Bayesian$Results$Results$mle$Z %>% head()
Latent.Bayesian$Results$Results$pmode$Z %>% head()
plot(Latent.Bayesian$Results$Results)
edges <- igraph::graph.adjacency(as.matrix(round(Social_Matrix[1:25,1:25]/10,0)), mode = "undirected", diag = F, weighted = T)
igraph::tkplot(edges, vertex.size = 30, edge.width = 1,vertex.label = colnames(Social_Matrix)[1:25],
               layout = Latent.Bayesian$Results$Results$mle$Z, vertex.color = I("Light Green"),
               vertex.label.cex = 1, vertex.label.color = "Black")

