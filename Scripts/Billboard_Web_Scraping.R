######################################################################
# Início : 03/04/2019                                                #
# Última modificação : 05/01/20                                      #
# Autor : Pedro Henrique Sodré Puntel                                #
# Email : pedro.puntel@gmail.com                                     #
# Instituição : Escola Nacional de Ciencias Estatísticas - ENCE IBGE #
# Disciplina : Projeto de Iniciação Científica - PIBIC CNPq 2018     #
# Tema : Web Scraping - https://www.billboard.com/charts/artist-100  #
######################################################################

##############
## Descrição :
############## 
# Neste scrpit, em um primeiro momento, será feito o web scraping com o site
# da Billboard (mais especificamente, da Artist 100 Chart) e posteriormente a
# montagem da sócio-matriz com os artistas associada.
#
# A chart mencionada é um ranqueamento dos Top 100 artistas da semana feito
# com base em diversos critérios (não exatamente claros) pela Billboard.
# Apesar da periodicidade semanal, existe observa-se pouca variação entre
# ranqueamentos produzidos entre semanas próximas, fato este que nos levou
# a considerar para fins de extração somente a última semana de cada mês,
# na esperaçnca de que esta contemple, em sua maior parte, toda variabilidade
# do mesmo.
#
# O procediemnto de Web Scraping para este site aproveita a boa estruturação
# e consistência da página https://www.billboard.com/charts/artist-100, no
# sentido de recriar por manipulação de strings as URLs das charts de acordo
# com o padrão observado e, posteriormente, acessar cada uma desta dentro de
# uma estrutura de repetição. De fato, o Web Scraping aqui perde um pouco do
# seu "brilho" por não fazer o uso de ferramentas mais avançadas como o Selenium,
# mas certamente ganha em perfomance.
#
# . Data da primeira chart produzida : 19/07/2014
# . Qualquer data selecionada anterior à 19/07/2014 é automaticamente redirecionada
# . Por exemplo, tente acessar: https://www.billboard.com/charts/artist-100/2000-01-01

##########
## Setup :
##########
# Pacotes utilizados
library("dplyr")
library("xml2")
library("rvest")
library("stringr")
library("lubridate")
library("compiler")

#################
## Web Scraping :
################# 
# Rotina que automatiza a coleta das informações e a montagem da sócio-matriz 
billboard_web_scraping = function() {
  
  # Estrutura básica da url Top 100 Artists Chart da Billboard
  chart_url <- "https://www.billboard.com/charts/artist-100"
  
  # Extraí data da chart mais recente
  beggining_date = read_html(chart_url) %>%
    html_node('#main > div.chart-detail-header > div.container.container--no-background.chart-detail-header__chart-info > div > span.dropdown.chart-detail-header__date-selector > button') %>%
    html_text(trim = TRUE) %>%
    str_remove_all(pattern = ',') %>% 
    as.Date("%B%d%Y")
  
  # Data da primeira chart produzida
  end_date = "2014-07-19" %>% as.Date()
  
  # Sequência de datas associadas as charts que iremos puxar (retrocede mensalmente no tempo)
  months_to_pull = seq(from = beggining_date, to = end_date, by = "-1 month")
  
  # Recriando as URLs
  charts_urls = paste0("https://www.billboard.com/charts/artist-100/", months_to_pull)
  
  # Objeto tipo lista que guardará os ranqueamentos das n-ésimas charts
  charts_list = list()
  
  # Navega e coleta os ranqueamentos
  for (i in 1:length(charts_urls)) {
    
    charts_list[[i]] = read_html(charts_urls[i]) %>%
      html_nodes(".chart-number-one__title , .chart-list-item__title") %>%
      html_text(trim = TRUE)
    
    # Necessário para evitar o erro HTTP 429
    Sys.sleep(3)
    
  }
  
  # Seleciona somente o Top 20 de cada chart eliminando a redund?ncia do primeiro artista
  charts_list = lapply(charts_list, function(i) {i = i[2:21]})
  
  # Vetor com todos os artistas que já ocuparam as charts (sem repetição)
  artists_unique = unlist(charts_list) %>% unique()
  
  # Todas as combinações de pares de artistas que deverão ser verificados para montagem da sócio-matriz
  artists_pairs = expand.grid(artists_unique, artists_unique)
  artists_pairs = with(artists_pairs, split(artists_pairs, as.factor(artists_pairs$Var2)))
  
  # Construção da sócio-matriz
  social_matrix = matrix(0, ncol = length(artists_unique), nrow = length(artists_unique))
  rownames(social_matrix) = artists_unique
  colnames(social_matrix) = artists_unique
  
  # Preenchimento da sócio-matriz
  for (i in 1:nrow(social_matrix)) {
    for (j in 1:ncol(social_matrix)) {
      
      counter = 0
      current_pair = c(artists_unique[i], as.character(artists_pairs[[i]][j,1]))
      
      counter = lapply(charts_list, function(k) {
        chart_artists  = as.character(k)
        flag = any(current_pair[1] == chart_artists) & any(current_pair[2] == chart_artists)
        if (isTRUE(flag)) {counter = counter + 1}
      })
      
      social_matrix[i,j] = sum(unlist(counter))
      counter = NULL
      
    }
  }
  
  # Retorna as informações
  out_list = list("Social_Matrix" = social_matrix, "Charts_List" = charts_list, "Unique_Artists" = artists_unique)
  
}
billboard_web_scraping = cmpfun(billboard_web_scraping)

#######################
## Obtenção dos dados :
#######################
# Inicia o processo
result = billboard_web_scraping()


