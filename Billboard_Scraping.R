######################################################################
# Data :  15/12/2018                                                 #
# Autor : Pedro Henrique Sodré Puntel                                #
# Instituição : Escola Nacional de Ciencias Estatísticas - ENCE IBGE #
# Disciplina : Projeto de Iniciação Científica - PIBIC CNPq 2018     #
# Tema : Top 100 Artists chart scraping - Billboard Web Site         #
#######################################################################

# Encoding default deste script é WINDOWS-1252
#
# Neste scrpit, faremos o scraping com o site da Billboard. Em um primeiro momento,
# puxaremos os dados referentes à chart Top 100 Artists e posteriormente, montaremos
# a sociomatiz associada.
#
# A chart mencionada é criada semanalmente pela BillBoard, porém, optou-se por scrapear
# somente uma única chart por mês dentro do intervalo de tempo considerado. A razão por
# trás disto é que pouca variação na ordenação dos artistas é observada de semana para
# semana.
# 
# Data máxima suportada pelo site da Billboad : 19/07/2014
# Qualquer data abaixo de 19/07/2014 é automaticamente redirecionada para 19/07/2014.
# Prova real : >> https://www.billboard.com/charts/artist-100/2000-01-01
#
# Em um primeiro momento, pensou-se que seria mais fácil scrapear os botões 'Last Week' 
# 'Next Week' da página. Porém, estes são na verdade elementos href, o que fica custoso
# computacionalmente via Selenium. Assim, optou-se pela 'força bruta', que consiste em 
# recriar as URLs das charts de acordo com o padrão observado e assim, acessar cada uma 
# destas para fazer a coleta. De fato, o scrap perde um pouco do seu "brilho" por não 
# fazer o uso de um ferramenta mais avançada como o Selenium, mas ganha em perfomance.

# Pacotes utilizados
library(usethis)
library(dplyr)
library(xml2)
library(rvest)
library(stringr)
library(lubridate)
library(foreach)
library(parallel)

# GitHub Setup :
use_git_config(scope = "user", user.name = "PedroPuntel", user.email = "pedro.puntel@gmail.com")

# Performance Setup :
cores_cluster = makeCluster(detectCores(), type = "PSOCK")
registerDoSEQ()

# Estrutura básica da Top 100 Artists Chart da Billboard
Chart.URL = "https://www.billboard.com/charts/artist-100"

# Extraí a data da chart mais recente e constrói então a URL que apota para tal
Beggining.Date = read_html(Chart.URL) %>%
  html_node('#main > div.chart-detail-header > div.container.container--no-background.chart-detail-header__chart-info > div > span.dropdown.chart-detail-header__date-selector > button') %>%
  html_text(trim = TRUE) %>%
  str_remove_all(pattern = ',') %>% as.Date("%B%d%Y")

# Primeira Top 100 Artists Chart da Billboard : 26/07/2014
End.Date = "2014-07-19" %>% as.Date()

# Constrói uma sequência de datas associadas as charts que iremos scrapear (retrocede 4 semanas no tempo)
Months.To_Pull = seq(from = Beggining.Date, to = End.Date, by = "-4 week")

# Recriando as URLs
Charts.URLs = paste0("https://www.billboard.com/charts/artist-100/", Months.To_Pull)

# Objeto tipo lista que guardará os artistas das n-ésimas charts
Charts.List = list()

# Processo de navegação e scrap das páginas (com paralelização)
foreach(i = 1:length(Charts.URLs)) %dopar% {
  print(paste0("Scrapeando página n° ", i))
  Charts.List[[i]] = read_html(Charts.URLs[i]) %>%
    html_nodes(".chart-number-one__title , .chart-list-item__title") %>%
    html_text(trim = TRUE)
}
stopCluster(cores_cluster)

# Selecionando somente os Top 20 artistas de cada chart eliminando a redundância do primeiro artista
Charts.List <- lapply(Charts.List, function(i) {i = i[2:21]})
Charts.List[[1]]

# Vetor com todos os artistas que já apareceram nas charts (sem repetição)
Artists.Unique = unlist(Charts.List) %>% unique()

# Combinações de pares de artistas que deverão ser verificados para montagem da sociomatriz
Artists.Pairs = expand.grid(Artists.Unique,Artists.Unique)
Artists.Pairs = with(Artists.Pairs, split(Artists.Pairs, as.factor(Artists.Pairs$Var2)))

# Constrói a sociomatriz
Social_Matrix = matrix(0, ncol = length(Artists.Unique), nrow = length(Artists.Unique))
rownames(Social_Matrix) = Artists.Unique
colnames(Social_Matrix) = Artists.Unique
Social_Matrix[1:5,1:5]

# Preenchimento da sociomatriz
for (i in 1:nrow(Social_Matrix)) {
  for (j in 1:ncol(Social_Matrix)) {
    counter = 0
    current_pair = c(Artists.Unique[i],as.character(Artists.Pairs[[i]][j,1]))
    counter = lapply(Charts.List, function(k) {
      chart_artists  = as.character(k)
      flag = any(current_pair[1] == chart_artists) & any(current_pair[2] == chart_artists)
      if (isTRUE(flag)) {counter = counter + 1}
    })
    Social_Matrix[i,j] = sum(unlist(counter))
    counter = NULL
  }
}
Social_Matrix[1:5,1:4]

# Exporta a sociomatriz
setwd("C:\\Users\\pedro\\Desktop\\R\\PIBIC 2018\\Analyses\\BillBoard\\Database")
rio::export(Social_Matrix, "15062019_BillboardTop20_SocialMatrix.csv", format = ",")

# Limpando Memória
rm(Artists.Pairs, Charts.List, Social_Matrix, Artists.Unique, Beggining.Date, Chart.URL,
   Charts.URLs, counter, current_pair, End.Date, i, j, Months.To_Pull, cores_cluster)
