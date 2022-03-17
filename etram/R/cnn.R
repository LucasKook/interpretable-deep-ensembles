
#' CNN for UTKFace data
#' @export
cnn_utkface <- function(output_shape = NULL, mbl = FALSE, ll_activation = "linear",
                        dropout_rate = 0.3, ll_bias = FALSE, input_shape = dim(im)[2:4], ...) {

  m <- keras_model_sequential() %>%

    layer_conv_2d(input_shape = input_shape,
                  filters = 16, kernel_size = c(3, 3), padding = "same", activation = "relu") %>%
    layer_dropout(rate = dropout_rate) %>%
    layer_conv_2d(filters = 32, kernel_size = c(3, 3), padding = "same", activation = "relu") %>%
    layer_dropout(rate = dropout_rate) %>%
    layer_max_pooling_2d(pool_size = c(2, 2)) %>%

    layer_conv_2d(filters = 32, kernel_size = c(3, 3), padding = "same", activation = "relu") %>%
    layer_dropout(rate = dropout_rate) %>%
    layer_conv_2d(filters = 32, kernel_size = c(3, 3), padding = "same", activation = "relu") %>%
    layer_dropout(rate = dropout_rate) %>%
    layer_max_pooling_2d(pool_size = c(2, 2)) %>%

    layer_conv_2d(filters = 32, kernel_size = c(3, 3), padding = "same", activation = "relu") %>%
    layer_dropout(rate = dropout_rate) %>%
    layer_conv_2d(filters = 32, kernel_size = c(3, 3), padding = "same", activation = "relu") %>%
    layer_dropout(rate = dropout_rate) %>%
    layer_max_pooling_2d(pool_size = c(2, 2)) %>%

    layer_flatten() %>%
    layer_dense(units = 500, activation = "relu") %>%
    layer_dropout(rate = dropout_rate) %>%
    layer_dense(units = 50, activation = "relu") %>%
    layer_dropout(rate = dropout_rate) %>%

    layer_dense(output_shape, activation = ll_activation, use_bias = ll_bias, ... = ...) %>%
    {if (mbl){
      layer_trafo_intercept()(.)
    } else .
    }
  return(m)
}

#' CNN for MNIST data
#' @export
cnn_mnist <- function(output_shape = NULL, mbl = FALSE, ll_activation = "linear",
                      ll_bias = FALSE, input_shape = dim(im)[2:4], ...) {

  m <- keras_model_sequential() %>%
    layer_conv_2d(input_shape = input_shape,
                  filters = 32, kernel_size = c(3, 3), activation = "relu") %>%
    layer_max_pooling_2d(pool_size = c(2, 2)) %>%
    layer_conv_2d(filters = 64, kernel_size = c(3, 3), activation = "relu") %>%
    layer_conv_2d(filters = 64, kernel_size = c(3, 3), activation = "relu") %>%
    layer_max_pooling_2d(pool_size = c(2, 2)) %>%
    layer_flatten() %>%
    layer_dense(100, activation = "relu") %>%

    layer_dense(output_shape, activation = ll_activation, use_bias = ll_bias, ... = ...) %>%
    {if (mbl){
      layer_trafo_intercept()(.)
    } else .
    }
  return(m)
}

#' CNN for stroke data
#' @export
cnn_stroke <- function(output_shape = NULL, mbl = FALSE, ll_activation = "linear",
                       dropout_rate = 0.3, ll_bias = FALSE, input_shape = dim(im)[2:5], ...) {

  m <- keras_model_sequential() %>%
    layer_conv_3d(input_shape = input_shape, filters = 32, # 32
                  kernel_size = c(3, 3, 3), padding = "same",
                  activation = "relu") %>%
    layer_max_pooling_3d(pool_size = c(2, 2, 2)) %>%
    layer_conv_3d(filters = 32,
                  kernel_size = c(3, 3, 3), padding = "same", activation = "relu") %>%
    layer_max_pooling_3d(pool_size = c(2, 2, 2)) %>%
    layer_conv_3d(filters = 64, kernel_size = c(3, 3, 3), padding = "same",
                  activation = "relu") %>%
    layer_max_pooling_3d(pool_size = c(2, 2, 2)) %>%
    layer_conv_3d(filters = 64, kernel_size = c(3, 3, 3), padding = "same",
                  activation = "relu") %>%
    layer_max_pooling_3d(pool_size = c(2, 2, 2)) %>%

    layer_flatten() %>%
    layer_dense(units = 128, activation = "relu") %>% # 128
    layer_dropout(rate = dropout_rate) %>%
    layer_dense(units = 128, activation = "relu") %>% # 128
    layer_dropout(rate = dropout_rate) %>%
    layer_dense(output_shape, activation = ll_activation, use_bias = ll_bias, ... = ...) %>%
    {if (mbl){
      layer_trafo_intercept()(.)
    } else .
    }
  return(m)
}

#' CNN for melanoma data
#' @export
cnn_melanoma <- function(output_shape = NULL, mbl = FALSE, ll_activation = "linear",
                         dropout_rate = 0.3, ll_bias = FALSE, input_shape = dim(im)[2:5], ...) {

  m <- keras_model_sequential() %>%
    layer_conv_2d(filters = 16, kernel_size = c(3, 3), activation = "relu",
                  input_shape = input_shape, padding = "same") %>%
    layer_dropout(dropout_rate) %>%
    layer_conv_2d(filters = 16, kernel_size = c(3, 3), activation = "relu",
                  padding = "same") %>%
    layer_dropout(dropout_rate) %>%
    layer_max_pooling_2d(pool_size = c(2, 2)) %>%
    layer_conv_2d(filters = 32, kernel_size = c(3, 3), activation = "relu",
                  padding = "same") %>%
    layer_max_pooling_2d(pool_size = c(2, 2)) %>%
    layer_dropout(dropout_rate) %>%
    layer_conv_2d(filters = 32, kernel_size = c(3, 3), activation = "relu",
                  padding = "same") %>%
    layer_dropout(dropout_rate) %>%
    layer_max_pooling_2d(pool_size = c(2, 2)) %>%

    layer_flatten() %>%
    layer_dense(units = 256, activation = "relu") %>%
    layer_dropout(dropout_rate) %>%
    layer_dense(units = 64, activation = "relu") %>%
    layer_dropout(dropout_rate) %>%
    layer_dense(output_shape, activation = ll_activation, use_bias = ll_bias, ... = ...) %>%
    {if (mbl){
      layer_trafo_intercept()(.)
    } else .
    }
  return(m)
}

