############################################################################
# Data :  15/03/2018                                                       #
# Autor : Pedro Henrique Sodré Puntel                                      #
# Instituição : Escola Nacional de Ciencias Estatísticas - ENCE IBGE       #
# Disciplina : Projeto de Iniciação Científica - PIBIC CNPq 2018/2019      #
# Professor : Gustavo Ferreira                                             #
# Tema : Billboard Top 100 Artists - Social Network Analysis               #
############################################################################

##########
## Notas :
########## 
## Default file encoding is WINDOWS-1252

################
## Referências :
################ 
## . Network Science
##   > https://www.sci.unich.it/~francesc/teaching/network/
##
## . Introduction to Social Network Methods 
##   > http://faculty.ucr.edu/~hanneman/nettext/
##
## . Social Network Analysis : A methodological introduction 
##   > http://courses.washington.edu/ir2010/readings/butts.pdf
##
## . Medidas de Centralidade em Grafos 
##   > http://objdig.ufrj.br/60/teses/coppe_m/LeandroQuintanilhaDeFreitas.pdf
## 
## . Comunity Analysis Algorithms in R using the igraph package
##   > https://stackoverflow.com/questions/9471906/what-are-the-differences-between-community-detection-algorithms-in-igraph#
##
## . R Tutorial: How to identify communities of items in networks :
##   > https://psych-networks.com/r-tutorial-identify-communities-items-networks/
## 
## . Introduction to R and Network Analysis :
##   > http://kateto.net/ruworkshop
##
## . Polnet Network Visualization with R :
##   > http://kateto.net/network-visualization

########################
## Pacotes necessários :
######################## 
library(dplyr)
library(igraph)
library(compiler)
library(data.table)
options(scipen = 999, digits = 4)

##########################
## Importa a Sociomatriz :
########################## 
f.path <- "C:\\Users\\pedro\\Desktop\\R\\PIBIC 2018\\Analyses\\BillBoard\\Database\\29052019_BillboardTop20_SocialMatrix.csv"
Social_Matrix <- rio::import(file = f.path) %>% as.matrix()
rownames(Social_Matrix) <- colnames(Social_Matrix)
rm(f.path)

## Redução da sociomatriz
Social_Matrix <- Social_Matrix[order(rowSums(Social_Matrix), decreasing = T),order(colSums(Social_Matrix), decreasing = T)]
Social_Matrix <- Social_Matrix[1:100,1:100]

####################
## Análise de Rede :
#################### 
## Grafo que modela a rede
Billboard.Graph <- graph_from_adjacency_matrix(Social_Matrix, mode = "undirected", diag = T, weighted = T)
V(Billboard.Graph)$name <- colnames(Social_Matrix)
is.weighted(Billboard.Graph)

## Função que computa características globais acerca do Grafo
Compute.Graph_Level_Indices <- function(g, use.loops, has.directions, has.weights) {
  
  ## Sub-rotina que fornece características básicas do grafo
  Graph.Stats <- function(g) {
    
    ## Número de arestas
    G.Edge_Count <- ecount(g)
    
    ## Número de vértices
    G.Vertex_Count <- vcount(g)
    
    ## Densidade do Grafo
    ## [!] Quão coeso/esparso é o grafo
    G.Density <- graph.density(g, loops = use.loops)
    
    ## Caminho Geodésico médio
    ## [!] Em média, quantas conexões cada vértice possui
    G.Avg_Path_Length <- mean_distance(g, directed = has.directions, unconnected = isFALSE(is.connected(g, mode = "strong")))
    
    ## Diâmetro do Grafo
    ## [!] Maior caminho geodésico
    G.Diameter <- farthest_vertices(g, directed = has.directions, unconnected = isFALSE(is.connected(g, mode = "strong")))
      diameter(g, directed = has.directions, unconnected = is.connected(g, mode = "strong"))
    
    ## Transitividade do Grafo
    ## [!] Quão interligadas são as conexões entre os vértices do grafo.
    ## [!] Na literatura, para grafos não-direcionados, a transitividade se dá pela razão entre
    ## o número de loops de tamanho 3 (x -- > y --> z -- > x) e o total de arestas do grafo
    if (isTRUE(has.directions)) {
      G.Transitivity <- transitivity(g, type = "directed", isolates = NaN)
    } else {
      G.Transitivity <- transitivity(g, type = "undirected", isolates = NaN)
    }
    
    ## Pareamento Assortativo/Dissortativo
    ## [!] Um fenômeno comum em redes sociais é a tendência natural dos atores em associar-se
    ## com aqueles semelhantes a si mesmos de alguma forma. Tal tendência é chamada de 
    ## homofilia ou mistura combinatória. É possível ainda, porém de forma mais rara, também
    ## encontramos uma mistura desassortativa, ou seja, a tendência dos se associarem com
    ## outros diferentes de si.
    ##  O pareamento assortativo de atores em uma rede pode ser identificado por meio de 
    ## características enumerativas ou escalares. Uma característica enumerativa tem um
    ## conjunto finito de valores possíveis. Exemplos são gênero, raça e nacionalidade.
    ## Dada uma característica enumerativa, cada nó da rede pode ser atribuído a um tipo de
    ## acordo com o valor da característica para o nó. Assim, dizemos que uma rede é então
    ## assortativa se uma fração significativa das arestas ligam entre vértices do mesmo tipo.
    ##  Matematicamente, uma medidade de assortatividade é definida pelo número de arestas
    ## que ocorrem entre vértices do mesmo tipo diminúidas daquelas as quais esperaríamos 
    ## encontrar em uma disposição aleatória dos vértices na rede. Indo além, pode-se provar 
    ## que a assortatividade de um vértice é dada pela sua covariância com os demais pares e
    ## que ao ser 'normalizada' coincide com o coeficiente de correalçaõ de Pearson.
    ## [!] Note que a definção de 'tipo' aqui se confunde com as medidas de centralidade.
    ## [!] https://stackoverflow.com/questions/31420117/r-and-igraph-help-assortativity-coefficient-with-weighted-edges-remaining-weig
    ## [!] https://arxiv.org/pdf/physics/0607134.pdf
    G.Assortativity <- assortativity_degree(g, directed = has.directions)
    
    ## Retorna as informações
    list("Edge_Count" = G.Edge_Count,
         "Vertex_Count" = G.Vertex_Count,
         "Density" = G.Density,
         "Avg_Geodesic" = G.Avg_Path_Length,
         "Diameter" = G.Diameter, 
         "Transitivity" = G.Transitivity,
         "Assortativity" = G.Assortativity) %>% return()
  }
  
  ## Sub-rotina que computa diversos índices de centralização do grafo
  Centralization.Index <- function(g) {
    
    ## [!] "Quão centralizado é o grafo em torno dos vértices com maior centralidade de (...)"
    
    ## Centralização por Grau
    G.Degree_Centralization <- centr_degree(g, mode = "all", loops = use.loops, normalized = T)
    
    ## Centralização por centralidade de proximidade
    G.Closeness_Centralization <- centr_clo(g, mode = "all", normalized = T)
    
    ## Centralização por centralidade de autovetor
    G.EigenVec_Centralization <- centr_eigen(g, directed = has.directions, normalized = T)
    
    ## Centralização por centralidade de intermediação
    G.Betweenness_Centralization <- centr_betw(g, directed = has.directions, normalized = T)
    
    ## Retorna as informações
    list("By_Degree" = G.Degree_Centralization,
         "By_Closeness" = G.Closeness_Centralization,
         "By_EigenVec" = G.EigenVec_Centralization,
         "By_Betweenness" = G.Betweenness_Centralization) %>% return()
  }
  
  ## Acessa as sub-rotinas e retorna os resultados
  list("Graph_Stats" = Graph.Stats(g), "Centralization_Index" = Centralization.Index(g)) %>% return()
  
}
Compute.Graph_Level_Indices <- cmpfun(Compute.Graph_Level_Indices)
Billboard.GLIs <- Compute.Graph_Level_Indices(Billboard.Graph, T, F, T)

## Análise das características globais do Grafo
##
## . 100 vértices
## . 2508 arestas
##
## . Densidade : 0.4966/1
##   --> Razoável. Do total de conexões que são teoricamente possíveis de existir entre
##       os atores da rede, ~50% realmente ocorrem. No nosso contexto, isso significa
##       que ~50% dos artistas mais populares da Billboard ocupam o Top 20 de forma
##       "eficiente".
##
##   --> Este resultado vai de encontro com os parâmetros estimados pelo Modelo de 
##       Abertura de conexões.
##
## . Caminho Geodésico Médio : 1.514 
##   --> Baixo. Em média, os vértices da rede estabelecem 1.5 conexões com os demais.
##       No nosso contexto, isso significa que os artistas mais populares da Billboard
##       compartilham, em média, 1.5 meses no Top 20 com os demais outros populares.
##
##   --> Uma outra interpretação possível para esta medida seria que dado o seu valor
##       baixo, existem muitos artistas dentre os populares que ocuparam o Top 20 de
##       forma esporádica (provavelmente por conta de uma música que bombou durante
##       um curto período de tempo). Isto pode ser verificado se inspecionarmos a 
##       diagonal da sociomatriz, por exemplo.
##
##   --> Vale ressaltar que este resultado é comum em redes socias no geral, devido
##       ao "Small World effect", conforme destacado em :
##           https://www.sci.unich.it/~francesc/teaching/network/geodesic.html
##
## . Diâmetro : 4 
##   --> Valor númerico baixo, porém expressivo no nosso contexto. O valor "4"
##       significa que a distância valorada entre os vértices mais afastados 
##       da rede é de 4 unidades. Ou seja, os artistas mais afastados da rede
##       compartilharam o Top 20 da Billboard 4 meses, algo que é no mínimo
##       interessante.
##
##   --> Porém, relacionando o diâmetro da rede com o caminho geodésico médio, o
##       o valor obtido é justificável uma vez que os artisas mais populares da 
##       Billboard compartilham uns com os outros, em média, 1.5 meses no Top 20.
##       
##   --> Assim, não é tão anormal pensar que os artisas mais afastados, Meghan Trainor
##       e WALK THE MOON compartilharam o Top 20 da Billboard por um período de 4 meses.
##
## . Transitividade : 0.6586/1
##   --> Razoável (se não alta). Significa dizer que ~65% das conexões entre os
##       vértices da rede são transitivas. No nosso contexto, isso equivale a dizer
##       que é possível "acessar" todos os artistas que mais ocuparam o Top 20 da
##       Billboard a partir de uma pequena parcela destes.
##
## . Assortatividade (por Centralidade de Grau) : -0.1269 de [-1,1]
##   --> Baixa e negativa. Em geral, é falsa a premissa de que vértices com alta
##       centralidade de Grau (ou seja, com muitos vizinhos) tendem a "assorciar-se"
##       uns com os outros. No nosso contexto, tal medida sugere que os artistas que
##       mais compartilharam o Top 20 com os demais outros não tendem a ficarem próximos
##       na rede, ou seja, parece que não existe um "grupinho dos mais populares dentre
##       os populares". 
## 
## . Centralização por Grau : 0.5115/1
##   --> Razoável. Espera-se uma grafo razoavelmente centralizado em torno dos vértices com
##       alta centralidade de grau, ou seja, em torno dos artistas que mais compartilharam
##       o Top 20 da Billboard com os demais.
##
## . Centralização por Proximidade : 0.6615/1
##   --> Razoável (se não alta). Espera-se um grafo centralizado em torno dos vértices com
##       alto valor de centralidade de proximidade, ou seja, em torno dos vértices cuja
##       distância geodésica dos demais outros é a baixa. No nosso contexto, isso equivale
##       a dizer que a rede composta pelos artistas mais populares da Billboard é razoavel-
##       mente centralizada em torno dos mais populares, dentre os populares. 
##
##   --> Naturalmente, como era de se esperar, essa medida vai de encontro com a anterior.
##
## . Centralização por centralidade de Autovetor : 0.4481/1
##   --> Mais uma vez, razoável. Espera-se um grafo em geral centralizado em torno dos
##       vértices com alta centralidade de autovetor, ou seja, em torno dos artistas que,
##       mais uma vez, compartilham o Top 20 com muitos outros.
##
##   --> Como esta já é a 3° medida que aponta para a mesma interpretação, é razoável dizer
##       que os "artistas mais importantes da Billboard" certamente ocuparão o centro do
##       grafo. Não obstante, essa premissa é confirmada se analisarmos o grafo referente
##       ao Modelo de Distância Latentes.
##  
## . Centralização por centralidade de Intermediação : 0.03254/1
##   --> Surpreendentemente baixa. Sugere totalmente o contrário das medidas anteriores,
##       não sei interpretar... Um palpite inicial que explicaria o seu baixo valor é
##       que a rede é  fortemente subdividida em grupos de artistas, mas não sei se é o
##       caso. Uma análise de comunidades seria necessária pra validar este palpite.

## Função que computa diversas medidas associadas a cada vértice do grafo
Compute.Node_Level_Indices <- function(g, use.loops, has.directions, has.weights) {
  
  ## Centralidade de Grau dos vértices do grafo
  ## . Número de vizinhos/vértices adjacentes de cada vértice do grafo
  Node.Deg <- degree(g, loops = use.loops, mode = "all", normalized = T)
  
  ## Centralidade do autovetor dos vértices do grafo
  ## . Nem todos os vértices de uma rede são equivalentes, alguns são mais importantes
  ## do que outros. Sendo assim, um vértice que estabelece algum tipo de relação com 
  ## um outro vértice "importante" da rede merece atenção se comparada com as demais.
  ## . "A node receiving many links does not necessarily have a high eigenvector centrality
  ## (it might be that all linkers have low or null eigenvector centrality). Moreover, a
  ## node with high eigenvector centrality is not necessarily highly linked (the node might
  ## have few but important linkers)."
  Node.Egv <- eigen_centrality(g, directed = has.directions, scale = T)
  
  ## Centralidade de Page Rank dos vértices do grafo
  ## . Um problema com a medida de centralidade de autovetor é que se um vértice com alta 
  ## centralidade referencia muitos outros, então estes últimos por sua vez também terão
  ## um alta centralidade. Porém, é razoável pensar que vértices referenciados por outros
  ## altamente centrais são de menor importância: o ganho de centralidade do referenciado
  ## deve ser menosprezado se quem o referencia, referencia também vários outros. Assim,
  ## a cálculo da centralidade Page Rank leva em consieração o número de arestas indicentes
  ## ao vértice, a propensão daquele que referencia, e por fim o valor da sua centralidade.
  ## . Variados estudos sugerem o valor 0.85 para o parâmetro 'damping' do Page Rank.
  Node.Pgrk <- page_rank(g, algo = "power", directed = has.directions, damping = 0.85)
  
  ## Centralidade de Kleingberg dos vértices do grafo
  ## . É razoável assumir também que um vértice é importante também por referenciar outros
  ## ivértices importantes da rede. Em redes de co-autorias, por exemplo, artigos podem
  ## referenciar vários outros os quais julga serem 'autoridades' no assunto e por isso,
  ## podem ser vistos como importantes. Seguindo essa lógica, é natural dividirmos os atores
  ## de uma rede em duas categorias: os 'chefes', que são referenciados por muitos outros,
  ## bem como aqueles 'subordinados', que apontam para os 'chefes'.
  Node.Klbrg <- authority.score(g, scale = T)
  
  ## Centralidade de Proximidade dos vértices do grafo
  ## . Computa a distância média que um vértice possui dos demais, retornando valores menores
  ## para os vértices que estão separados por uma curta distância geodésica. Vértices com
  ## essa característica possuem um melhor acesso a informação dissipada na rede. Em redes
  ## sociais, atores com menores índices de centralidade de proximidade são mais propícios
  ## a influenciar os demais vértices.
  Node.Clo <- closeness(g, mode = "all", normalized = T)
  
  ## Centralidade de Entrelaçamento dos vértices do grafo
  ## . Quantifica o quão provável um dado vértice é de intermediar a ligação entre
  ## outros dois vértices. Nesse sentido, vértices com alto índice de entrelaçamento
  ## influenciam diretamente no fluxo de informação pela rede. Analagomente, estes são
  ## aqueles cuja remoção da rede irá pertubar a comunição inter-vértices.
  Node.Btw <- betweenness(g, directed = has.directions, normalized = T)
  
  ## Data Table com as medidas de centralidade de cada vértice da rede
  Degree <- Node.Deg %>% as.list() %>% unlist() %>% as.numeric()
  EigenVec <- Node.Egv$vector %>% as.list() %>% unlist() %>% as.numeric()
  PageRank <- Node.Pgrk$vector %>% as.list() %>% unlist() %>% as.numeric()
  Kleinberg <- Node.Klbrg$vector %>% as.list() %>% unlist() %>% as.numeric()
  Closeness <- Node.Clo %>% as.list() %>% unlist() %>% as.numeric()
  Betweenness <- Node.Btw %>% as.list() %>% unlist() %>% as.numeric()
  Node.Centrality_Table <- cbind(Degree, EigenVec, PageRank, Kleinberg, Closeness, Betweenness) %>% as.matrix()
  rownames(Node.Centrality_Table) <- V(g)$name
  
  ## Retorna a tabela
  list("Centrality_Tale" = Node.Centrality_Table) %>% return()
  
}
Compute.Node_Level_Indices <- cmpfun(Compute.Node_Level_Indices)
Billboard.NLIs <- Compute.Node_Level_Indices(Billboard.Graph, T, F, T)
tb <- Billboard.NLIs$Centrality_Tale

## Centralidade de Grau
## . Sugere uma clara categorização dos artistas
##   --> Consistentemente Populares : Drake, Taylor Swift, ..., Justin Bieber
##   --> Os que ainda dão pro gasto : Meghan Trainor, Rihana, ..., J. Cole
##   --> Aqueles que volta e meia bombam : Fetty Wap, Migos, ..., Gucci Mane
##   --> Aqueles que brotam do nada : Iggy Azaela, ..., Jonas Brothers
## . Distribuição aparentemente assimétrica à direita, sugere vai uma característica
## clássica da análise de redes sociais : "Muitos são os atores em uma rade com poucos
## vizinhos e poucos são aqueles com muitos vizinhos."
tb[,1] %>% as.numeric() %>% hist(main = "Distribuição - Centralidade de Grau",
                                 ylab = "Frequência", freq = T, breaks = "Sturges",
                                 xlab = "Centralidade de Grau", col = "Gray",
                                 border = "Blue")
tb[order(tb[,1], decreasing = T),1] %>% head(n = 10L)
tb[order(tb[,1], decreasing = T),1] %>% tail(n = 10L)

## Centralidade de Autovetor
## . Sugere uma ordenção levemente diferente dos artistas. Porém, ainda sim, convêm 
## confiarmos mais nos resultados dessa medida do que na centraldiade de Grau por conta
## do seu embasamento teórico "mais robusto" (ler descrição dentro da função).
## . Distribuição mais acentuadamente assimétrica do que a centralidade de Grau, reforça
## a ideia de que poucos são os artistas "extraordinários" na rede.
tb[,2] %>% as.numeric() %>% hist(main = "Distribuição - Centralidade de Autovetor",
                                 ylab = "Frequência", freq = T, breaks = "Sturges",
                                 xlab = "Centralidade de Autovetor", col = "Gray",
                                 border = "Blue")
tb[order(tb[,2], decreasing = T),2] %>% head(n = 10L)
tb[order(tb[,2], decreasing = T),2] %>% tail(n = 10L)

## Centralidade de Page Rank
## . Sugere uma ordenação dos artistas parecida com a de autovetor.
## . Destaque para a Beyonce, Sam Smith e Florida Geogia lina que subiram muito.
## . Interessante notar que intervalo de variação dos valores para a centralidade
## são bem diferentes daqueles observados nas demais medidas. Em um primeiro momento,
## não sei dizer porquê isso ocorre, mas é possível que a centralidade de Page Rank
## não seja capaz de diferenciar os artistas na rede dada a pequena variação de seus
## valores.
## . Distribuição mais parecida com aquela da centralidade de Grau. Isto é um tanto
## que inesperado, uma vez que a centralidade de Page Rank deriva diretamente da de
## Autovetor.
tb[,3] %>% as.numeric() %>% hist(main = "Distribuição - Centralidade de Page Rank",
                                 ylab = "Frequência", freq = T, breaks = "Sturges",
                                 xlab = "Centralidade de Page Rank", col = "Gray",
                                 border = "Blue")
tb[order(tb[,3], decreasing = T),3] %>% head(n = 10L)
tb[order(tb[,3], decreasing = T),3] %>% tail(n = 10L)

## Centralidade de Kleinberg
## . Ordenação dos artistas e distribuição parecida com a de Autovetor.
## . Maior variabiliade da medida (consegue diferenciar melhor os artistas)
## . Até o momento, esta é a medida de centralidade que melhor retrata a rede.
tb[,4] %>% as.numeric() %>% hist(main = "Distribuição - Centralidade de Kleinberg",
                                 ylab = "Frequência", freq = T, breaks = "Sturges",
                                 xlab = "Centralidade de Kleinberg", col = "Gray",
                                 border = "Blue")
tb[order(tb[,4], decreasing = T),4] %>% head(n = 10L)
tb[order(tb[,4], decreasing = T),4] %>% tail(n = 10L)

## Centralidade de Proximidade
## . Como era de se esperar, sugere uma ordenação totalmente distinta dos artistas
## na rede. Isso se deve porque o propósito da centralidade de proximidade é quantificar
## o quão próximo um vértice na rede é distante dos demais. 
## . Assim, uma explicação lógica (e correta) para tal discrepância é que os artistas
## Logic, Gucci Mane, etc... estão todos tão intensamente próximos uns dos outros que
## a sua distância aos demais vértices chega a ser "irrelevante" em termos de cálculo.
## Não obstante, observe que os nossos artistas superpopulares, Drake, Ed Sheeran, etc...
## são justamente aqueles com o menor valor desta medida, pois estão todos isolados dos
## demais.
## . Indiretamente, essa medida sugere uma forte presença de comunidades de artistas,
## pois somente um cluster densamente conectado seria capaz de anular a distância aos
## demais outros vértices da rede. Como veremos a frente, tal hipótese é facilemnte 
## confirmada na Análise de Comunidades e na inspeção do Grafo.
tb[,5] %>% as.numeric() %>% hist(main = "Distribuição - Centralidade de Proximidade",
                                ylab = "Frequência", freq = T, breaks = "Sturges",
                                xlab = "Centralidade de Proximidade", col = "Gray",
                                border = "Blue")
tb[order(tb[,5], decreasing = T),5] %>% head(n = 10L)
tb[order(tb[,5], decreasing = T),5] %>% tail(n = 10L)

## Centralidade de Intermediação
## . Mais uma vez, como era de se esperar, uma ordenação totalmente diferente dos
## aristas. Porém, de forma similar a centralidade de proximidade, a centralidade de
## intermediação tem um propósito muito específico : ela quantifica o quanto que um
## vértice da rede é responsável por conectar todos os demais uns com os outros.
## . A medida atribuí valores altos para artistas como Prince (não sei nem quem é),
## Kanye West, Bruno Mars, etc... pois estes artisas "intermediam fortemente o fluxo de
## informação" dentre a comunidade que pertencem, ou seja, remover estes artistas da rede
## é "destruir a ecleticidade dos artistas que ocuparam o Top 20 ao seu redor".
## . Note que artistas como Drake, Ed Sheeran, etc... são justamente aqueles com os
## menores valores desta medida.
tb[,6] %>% as.numeric() %>% hist(main = "Distribuição - Centralidade de Intermediação",
                                 ylab = "Frequência", freq = T, breaks = "Sturges",
                                 xlab = "Centralidade de Intermediação", col = "Gray",
                                 border = "Blue")
tb[order(tb[,6], decreasing = T),6] %>% head(n = 10L)
tb[order(tb[,6], decreasing = T),6] %>% tail(n = 10L)

###########################
## Análise de Comunidades :
########################### 

## Função que implementa diversas rotinas para detecção de comunidades
Compute.Graph_Communities <- function(g, use.loops, has.directions, has.weights) {
  
  ## Identificação de comunidades de atores em uma rede por maximização de modularidade
  ## . A modularidade de uma rede é uma medida que computa o número de arestas que ocorrem
  ## entre vértices que pertencem a uma mesma comunidade diminuída do número de arestas que
  ## esperaríamos encontrar em uma disposição aleatória dos vértices pela rede.
  ## . Tal medida assume valores positivos se existem mais arestas conectando vértices 'do
  ## mesmo tipo' do que o esperado e valores negativos caso contrário. Os algoritmos aqui 
  ## apresentados buscam um particionamento da rede em comunidades de vértices de forma a
  ## maximizar a modularidade.

  ## Método Fast Greedy
  ## . Método de abordagem hierárquica, porém que atua 'de baixo para cima'. Atua otimizando
  ## diretamente a função modularidade associada a rede de forma gananciosa. Inicialmente,
  ## cada nó pertence a uma comunidade separada, e as comunidades são mescladas iterativamente.
  ## Cada fusão é sempre localmente ótima (produz o maior aumento no valor da modularidade).
  ## O algoritmo termina quando não é mais possível aumentar a modularidade. O método é rápido,
  ## e geralmente tentado como uma primeira aproximação porque não possui parâmetros para
  ## ajustar. No entanto, sabe-se que este sofre de um limite de resolução, isto é, comunidades
  ## abaixo de um determinado limite de tamanho serão sempre mescladas com comunidades vizinhas.
  Fast_Greedy <- cluster_fast_greedy(g, merges = T, membership = T)
  V(g)$Fast_Greedy_Community_Attr <- Fast_Greedy$membership
  
  ## Método do Autovetor Líder
  ## . Método de abordagem hierárquica comum que otimiza a função de modularidade. Em cada etapa,
  ## o grafo é dividido em duas partes onde tal separação produz sempre um aumento de modularidade.
  ## Porém, aqui a divisão é determinada pela avaliação do autovetor associado ao maior autovalor da
  ## chamada matriz de modularidade. Devido aos cálculos do autovetor envolvidos, tal método pode não
  ## funcionar em gráficos degenerados. A vantagem deste se comparado aos demais é que para grafos
  ## não degenerados, atinge-se um score de modularidade maior do que o método fast greedy, embora
  ## pouco mais lento. Agrupa os vértices baseando-se nos sinais das entradas deste mesmo autovetor.
  ## Caso todas as entradas do autovetor sejam de mesmo sinal, significa que a rede possui 
  ## estruturação única.
  Egv_community <- cluster_leading_eigen(g, steps = 10, callback = NULL)
  V(g)$Egv_Community_Attr <- Egv_community$membership
  
  ## Método da Caminhanda Aleatória
  ## . Tenta encontrar comunidades em um grafo pelo método 'Random Walk'. A premissa deste
  ## método é que caminhadas curtas tendem a permanecer na mesma comunidade.
  Rdm_Walk <- cluster_walktrap(g, steps = 10, modularity = T, merges = T)
  V(g)$Rdm_Walk_Community_Attr <- Rdm_Walk$membership

  ## Retorna os grafos associados a cada algoritmo de clusterização
  g.layout <- layout_with_mds(g, dim = 2)
  
  n.cl <- as.factor(V(g)$Fast_Greedy_Community_Attr) %>% levels() %>% as.numeric() %>% max()
  g.clrs <- RColorBrewer::brewer.pal(n.cl, name = "Set1")
  plot.igraph(g, size = 25, vertex.shape = "none", main = "Comunity Analysis - Fast Greedy", 
              vertex.label = V(g)$name, vertex.label.font = 2, vertex.label.cex = 1,
              vertex.label.color = g.clrs[V(g)$Fast_Greedy_Community_Attr],
              edge.color = "White", edge.width = 0.75, edge.lty = "solid", layout = g.layout)
  
  n.cl <- as.factor(V(g)$Egv_Community_Attr) %>% levels() %>% as.numeric() %>% max()
  g.clrs <- RColorBrewer::brewer.pal(n.cl, name = "Set1")
  plot.igraph(g, size = 25, vertex.shape = "none", main = "Comunity Analysis - Leading Eigen Vector", 
              vertex.label = V(g)$name, vertex.label.font = 2, vertex.label.cex = 1,
              vertex.label.color = g.clrs[V(g)$Egv_Community_Attr], edge.color = "white",
              edge.width = 0.75, edge.lty = "solid", layout = g.layout)
  
  n.cl <- as.factor(V(g)$Rdm_Walk_Community_Attr) %>% levels() %>% as.numeric() %>% max()
  g.clrs <- RColorBrewer::brewer.pal(n.cl, name = "Set1")
  plot.igraph(g, size = 25, vertex.shape = "none", main = "Comunity Analysis - Random Walk", 
              vertex.label = V(g)$name, vertex.label.font = 2, vertex.label.cex = 1,
              vertex.label.color = g.clrs[V(g)$Rdm_Walk_Community_Attr], edge.color = "white",
              edge.width = 0.75, edge.lty = "solid", layout = g.layout)
  
  ## Retorna as informações dos algoritmos
  list("Fast_Greedy" = Fast_Greedy, "Lead_Eigen" = Egv_community, "Random_Walk" = Rdm_Walk) %>% return()
  
}
Compute.Graph_Communities <- cmpfun(Compute.Graph_Communities)
Billboard.Communities <- Compute.Graph_Communities(Billboard.Graph, F, F, T)

## Conclusões
## . Infelizmente, os algoritmos para detecção de comunidades retornam resultados
## contraditórios : alguns subdividem "de mais" a rede (que é o caso com o método
## Fast Greedy e Leading Eigenvector) enquanto outros subdividem menos a rede (como
## é o caso do Random Walk).
## . Não só o número de partições da rede que são diferentes entre os métodos, a própria
## partição em si não vai de encontro com os valores e as interpretações das medidas
## de centralidade. Por exemplo :
##
##  --> Porque Drake, Ed Sheraan, Shawn Mendes (ou seja, os mais top's) não
##      fazem parte de uma única comunidade ? Estes encontram-se fisicamente
##      próximos na rede, porém são subdivididos em grupos diferentes...
##
##  --> Não consigo pensar em uma razão, ou melhor, um fator que justifica o agrupamento
##      dos artistas... "Em um mesmo grupo, existem misturados artistas famosos bem como
##      artistas não tão famosos".
##
## . A confusão se extende mais ainda pois aparentemente, o Grafo não aparenta ser centra-
## lizado em torno dos vértices com os maiores índices de centraldiade de proximidade...
## Este é, no entanto, altamente centralizado em torno dos vértices com alto valor de 
## centralidade de Grau, Autovetor e Kleinberg. (Note como Drake, Ariana Grande estão
## de fato misturados no centro).

##################
## Visualizações :
################## 

## Layouts de disposição dos Vértices
##
## . "Random" : Disposição aleatória dos vértices, porém que respeita porém as suas
## respectivas centralidade de Grau.
##
## . "Eigen" : Dispõem os vértices com base nas estruturas espectrais da matriz
## de adjacência da rede.
##
## . "Fruchterman & Reingold / Kamada & Kawai" : Algoritmos que buscam simular uma
## interação de ordem física entre os vértices de uma rede visando determinar suas
## 'posições' ótimas. Aqui, vértices exercem uma força repulsiva uns com os outros
## e as arestas que os interligam atuam como a força contrária que os aproximam. A
## inspiração física se dá por analogia dos vértices como particulas eletricamente
## carregadas que se repelem e as arestas como molas que os atraem. Tais forças
## concomitam para um movimento de atração/repulsão que eventualmente convergirá 
## para um estado de equilíbro. Cada etapa do algoritmo é de tamanho pré-determinado
## (niter) e ao seu término, diz-se que a energia do sistema é minimizada.
##
## . MDS : Dispõe os vértices da no grafo de acordo com a semelhança no seu padrão
## de conexões.

## Parâmetros aceitos para os Layouts
## . random(dist = "normal")
## . mds(var = "geodist", dist = "euclidean")  
## . kamadakawai(niter = 1500)                 
## . fruchtermanreingold(niter = 1500)         
## . eigen(var = "symlower", evsel = "size")

## Grafos somente com os vértices 
par(mfcol = c(1,2))
net.lay <- layout_with_mds(Billboard.Graph, dim = 2)
V(Billboard.Graph)$name <- colnames(Social_Matrix)
plot.igraph(Billboard.Graph, main = "Billboard Network - Artists",
     vertex.shape = "none", edge.color = I("white"), edge.lty = "solid",
     vertex.label = V(Billboard.Graph)$name, vertex.label.cex = 0.85,
     layout = net.lay,
     vertex.label.color = I("Black"), vertex.label.font = 1,
     sub = "Multidimensional Scaling Layout")
V(Billboard.Graph)$name <- NA
plot.igraph(Billboard.Graph, main = "Billboard Network - Structure",
     vertex.shape = "circle", edge.color = I("White"),
     vertex.color = I("Dark Orange"), labelS = NA,
     layout = net.lay,
     sub = "Multidimensional Scaling Layout")

par(mfcol = c(1,2))
net.lay <- layout_with_fr(Billboard.Graph, dim = 2, niter = 2000)
V(Billboard.Graph)$name <- colnames(Social_Matrix)
p <- plot.igraph(Billboard.Graph, main = "Billboard Network - Artists",
            vertex.shape = "none", edge.color = I("white"), edge.lty = "solid",
            vertex.label = V(Billboard.Graph)$name, vertex.label.cex = 0.85,
            layout = net.lay,
            vertex.label.color = I("Black"), vertex.label.font = 1,
            sub = "Fruchterman & Reingold Layout")
V(Billboard.Graph)$name <- NA
q <- plot.igraph(Billboard.Graph, main = "Billboard Network - Structure",
            vertex.shape = "circle", edge.color = I("White"),
            vertex.color = I("Dark Orange"), labelS = NA,
            layout = net.lay,
            sub = "Fruchterman & Reingold Layout")

par(mfcol = c(1,2))
net.lay <- layout_with_kk(Billboard.Graph)
V(Billboard.Graph)$name <- colnames(Social_Matrix)
plot.igraph(Billboard.Graph, main = "Billboard Network - Artists",
            vertex.shape = "none", edge.color = I("white"), edge.lty = "solid",
            vertex.label = V(Billboard.Graph)$name, vertex.label.cex = 0.85,
            layout = net.lay,
            vertex.label.color = I("Black"), vertex.label.font = 1,
            sub = "Kamada & Kawai Layout")
V(Billboard.Graph)$name <- NA
plot.igraph(Billboard.Graph, main = "Billboard Network - Structure",
            vertex.shape = "circle", edge.color = I("White"),
            vertex.color = I("Dark Orange"), labelS = NA,
            layout = net.lay,
            sub = "Kamada & Kawai Layout")

## Grafos Normais
par(mfcol = c(1,1))
V(Billboard.Graph)$name <- colnames(Social_Matrix)
plot(Billboard.Graph, main = "Billboard Network",
     vertex.shape = "circle", vertex.color = I("Dark Orange"), vertex.size = 20,
     vertex.label = V(Billboard.Graph)$name, vertex.label.font = 1, vertex.label.cex = 0.9,
     vertex.label.color = I("Black"), edge.color = I("Grey"), edge.width = 1.25, edge.lty = 1,
     layout = layout_with_mds(Billboard.Graph, dim = 2))

plot(Billboard.Graph, main = "Billboard Network - Fruchterman & Reingold Layout",
     vertex.shape = "circle", vertex.color = I("Dark Orange"), vertex.size = 20,
     vertex.label = V(Billboard.Graph)$name, vertex.label.font = 2, vertex.label.cex = 0.9,
     vertex.label.color = I("Black"), edge.color = I("Grey"), edge.width = 1.25, edge.lty = 1,
     layout = layout_with_fr(Billboard.Graph, niter = 2000))

plot(Billboard.Graph, main = "Billboard Network - Kamada & Kawai Layout",
     vertex.shape = "circle", vertex.color = I("Dark Orange"), vertex.size = 20,
     vertex.label = V(Billboard.Graph)$name, vertex.label.font = 2, vertex.label.cex = 0.9,
     ertex.label.color = I("Black"), edge.color = I("Grey"), edge.width = 1.25, edge.lty = 1,
     layout = layout_with_kk(Billboard.Graph, niter = 2000))


