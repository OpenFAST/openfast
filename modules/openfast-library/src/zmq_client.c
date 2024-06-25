// ||  ----- ZeroMQ Client (re-adapted from ROSCO) ------ ||
// ||            Real-time interactor for FAST            ||
// || --------------------------------------------------- ||


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include </opt/homebrew/Cellar/zeromq/4.3.5_1/include/zmq.h>
#include </Users/lorenzoschena/Desktop/Workspace/openfast/cJSON/cJSON.h>

// Initializing publisher and requester as global variables

void *publisher = NULL;
void *requester = NULL;

void printString(const char *str) {
    printf("String: ");
    
    // Print each character of the string including the null terminator
    for (int i = 0; str[i] != '\0'; ++i) {
        printf("%c", str[i]);
    }
    printf("\\0");
}


// float *zmq_req_rep(const char *socket_address, const char *request) {
//     // add null-termination from Fortran strings inputs 

//     // printf("received request from fortran in C: %s \n\n", request);

//     size_t socket_len = strlen(socket_address);
//     size_t request_len = strlen(request);

//     int send_status_req;
//     int send_status_data;
//     // Allocate memory for null-terminated strings
//     // char *socket_copy = malloc(socket_len + 1);
//     // char *request_copy = malloc(request_len + 1);

//     // // Copy received strings and null-terminate them
//     // strncpy(socket_copy, socket_address, socket_len);
//     // socket_copy[socket_len +1] = '\0';

//     // strncpy(request_copy, request, request_len);
//     // request_copy[request_len +1] = '\0';

//     void *context = zmq_ctx_new();

//     if (context == NULL) {
//         perror("Error in opening ZMQ comm");
//         return NULL; // Return NULL to indicate failure
//     }

//     void *requester = zmq_socket(context, ZMQ_REQ);
//     if (requester == NULL) {
//         perror("Error in opening ZMQ comm");
//         zmq_ctx_destroy(context);
//         return NULL; // Return NULL to indicate failure
//     }
//     // printf("Connecting...");

//     int rc = zmq_connect(requester, socket_address);

//     if (rc != 0) {
//         perror("Error in connecting to specified address");
//         fprintf(stderr, "Address: %s\n", socket_address);
//         zmq_close(requester);
//         zmq_ctx_destroy(context);
//         return NULL; // Return NULL to indicate failure
//     }

//     // printf("C: Sending request %s to socket... \n", request);
//     send_status_req = zmq_send(requester, request, strlen(request), 0);

//     if (send_status_req >= 0) {
//         printf("C: Request sent successfully.\n");
//     } else {
//         printf("C: Error sending request: %s\n", zmq_strerror(errno));
//     }

//     float *received_value = (float *)malloc(sizeof(float));

//     int recv_size = zmq_recv(requester, (char *)received_value, sizeof(float), 0);
//     if (recv_size == sizeof(float)) {
//         printf("C: Received float value: %f\n", *received_value);
//         return received_value; // Return pointer to received float data

//     } else {
//         free(received_value); // Free the allocated memory in case of failure
//         printf("C: Error receiving float value\n");
//         zmq_close(requester);
//         zmq_ctx_destroy(context);

//         return NULL;
//     }

//     zmq_close(requester);
//     zmq_ctx_destroy(context);
//     free(socket_address);
//     free(request);

//     return NULL; // Return NULL to indicate failure
// }
// String manipulation before passing to ZMQ


char **split(const char *str, char delimiter, int *count) {
    int token_count = 0;
    int str_len = strlen(str);

    // Count the number of tokens
    for (int i = 0; i < str_len; i++) {
        if (str[i] == delimiter) {
            token_count++;
        }
    }
    token_count++;  // Add one for the last token

    char **tokens = malloc(token_count * sizeof(char *));
    if (tokens == NULL) {
        *count = 0;
        return NULL;
    }

    int token_index = 0;
    int token_start = 0;
    for (int i = 0; i <= str_len; i++) {
        if (str[i] == delimiter || str[i] == '\0') {
            int token_length = i - token_start;
            tokens[token_index] = malloc((token_length + 1) * sizeof(char));
            if (tokens[token_index] == NULL) {
                *count = token_index;
                for (int j = 0; j < token_index; j++) {
                    free(tokens[j]);
                }
                free(tokens);
                return NULL;
            }
            strncpy(tokens[token_index], str + token_start, token_length);
            tokens[token_index][token_length] = '\0';
            token_index++;
            token_start = i + 1;
        }
    }

    *count = token_count;
    return tokens;
}

char *createJSONString(const char *data, const char *names) {
    cJSON *root = cJSON_CreateObject();
    
    int data_count, names_count;
    char **data_tokens = split(data, ';', &data_count);
    char **names_tokens = split(names, ';', &names_count);

    if (root && data_tokens && names_tokens && data_count == names_count) {
        for (int i = 0; i < data_count; i++) {
            // Check for empty string in names before adding to JSON
            if (strcmp(names_tokens[i], "") != 0) {
                cJSON_AddItemToObject(root, names_tokens[i], cJSON_CreateNumber(atof(data_tokens[i])));
            }
            free(data_tokens[i]);
            free(names_tokens[i]);
        }
        free(data_tokens);
        free(names_tokens);

        char *jsonString = cJSON_Print(root);
        cJSON_Delete(root);
        return jsonString;
    } else {
        cJSON_Delete(root);
        if (data_tokens) {
            for (int i = 0; i < data_count; i++) {
                free(data_tokens[i]);
            }
            free(data_tokens);
        }
        if (names_tokens) {
            for (int i = 0; i < names_count; i++) {
                free(names_tokens[i]);
            }
            free(names_tokens);
        }
        return NULL;
    }
}

char *cJSON_keys_to_string(cJSON *json, const char *delimiter) {
    char *result = NULL;
    cJSON *child = json->child;
    size_t total_len = 0;

    // Calculate the total length needed for the keys string
    while (child != NULL) {
        total_len += strlen(child->string) + strlen(delimiter) + 1; // +1 for delimiter

        child = child->next;
    }

    // Allocate memory for the keys string
    result = (char *)malloc(total_len + 1); // +1 for null terminator
    if (result == NULL) {
        return NULL;
    }

    // Build the keys string separated by the delimiter
    child = json->child;
    result[0] = '\0'; // Initialize the string
    while (child != NULL) {
        strcat(result, child->string);
        if (child->next != NULL) {
            strcat(result, delimiter);
        }
        child = child->next;
    }

    return result;
}

char *cJSON_values_to_string(cJSON *json, const char *delimiter) {
    char *result = NULL;
    cJSON *child = json->child;
    size_t total_len = 0;

    // Calculate the total length needed for the values string
    while (child != NULL) {
        total_len += strlen(cJSON_Print(child)) + strlen(delimiter) + 1; // +1 for delimiter

        child = child->next;
    }

    // Allocate memory for the values string
    result = (char *)malloc(total_len + 1); // +1 for null terminator
    if (result == NULL) {
        return NULL;
    }

    // Build the values string separated by the delimiter
    child = json->child;
    result[0] = '\0'; // Initialize the string
    while (child != NULL) {
        strcat(result, cJSON_Print(child));
        if (child->next != NULL) {
            strcat(result, delimiter);
        }
        child = child->next;
    }

    return result;
}

float *process_received_data(const char *received_data) {
    // Parse received JSON data
    cJSON *root = cJSON_Parse(received_data);
    if (!root) {
        fprintf(stderr, "Error parsing JSON\n");
        return NULL; // Handle the error appropriately
    }

    // Extract values from JSON and create an array of floats
    cJSON *array = cJSON_GetObjectItem(root, "array_key"); // Replace "array_key" with your JSON key
    if (!array || !cJSON_IsArray(array)) {
        fprintf(stderr, "JSON array not found or invalid\n");
        cJSON_Delete(root);
        return NULL; // Handle the error appropriately
    }

    int num_elements = cJSON_GetArraySize(array);
    float *float_array = (float *)malloc(num_elements * sizeof(float));

    for (int i = 0; i < num_elements; ++i) {
        cJSON *item = cJSON_GetArrayItem(array, i);
        if (cJSON_IsNumber(item)) {
            float_array[i] = (float)item->valuedouble;
        } else {
            fprintf(stderr, "Non-numeric value found in JSON array\n");
            free(float_array);
            cJSON_Delete(root);
            return NULL; // Handle the error appropriately
        }
    }

    cJSON_Delete(root);
    return float_array;
}


// ------------------------------- Key ZMQ Helpers ------------------------- //
int zmq_init_pub(const char *req_address) {
    void *context = zmq_ctx_new();
    publisher = zmq_socket(context, ZMQ_PUB);

    int rc = zmq_connect(publisher, req_address);

    if (rc != 0) {
        printf("Error binding at %s : %s\n", req_address, zmq_strerror(errno));
        zmq_close(publisher);
        zmq_ctx_destroy(context);
        return -1; // Return error code if binding fails
    }

    printf("Established PUB - SUB connection at address %s \n", req_address);
    return 0; // Return success
}


int zmq_init_req_rep(const char *req_address) {
    void *context = zmq_ctx_new();
    requester = zmq_socket(context, ZMQ_REQ);

    int rc = zmq_connect(requester, req_address);

    if (rc != 0) {
        printf("Error binding at %s : %s\n", req_address, zmq_strerror(errno));
        zmq_close(requester);
        zmq_ctx_destroy(context);
        return -1; // Return error code if binding fails
    }

    printf("Established REQ - REP connection at address %s \n", req_address);
    return 0; // Return success
}

int zmq_broadcast(const char *data, const char *names) {
    if (publisher == NULL) {
        printf("Socket not initialized. Call zmq_initialize_publisher first.\n");
        return -1; // Return error if socket is not initialized
    }
    // printf("inside C broadcasting routine.... \n");

    size_t size = 0;

    // printf("Received names: %s\n", names);
    // printf("Received string: %s\n", data);

    char *jsonString = createJSONString(data, names);

    if (jsonString) {
        // printf("JSON String ready to be published: %s\n", jsonString);
        zmq_send(publisher, jsonString, strlen(jsonString), 0);
        free(jsonString); 
    } else {
        printf("Failed to create JSON string\n");
    }

    return 0; // Success
}

int count_semicolons(const char *request) {
    int semicolon_count = 0;
    for (int i = 0; request[i] != '\0'; i++) {
        if (request[i] == ';') {
            semicolon_count++;
        }
    }
    return semicolon_count;
}


float *zmq_req_rep(const char *socket_address, const char *request) {
    
    if (requester == NULL) {
        printf("Socket not initialized. Call zmq_initialize_requester first.\n");
        return NULL; // Return error if socket is not initialized
    }

    int req_count, send_status_req;
    int max_floats = count_semicolons(request); 

    // char **reqtoken = split(request, ';', &req_count);

    // printf("C: Sending request %s to socket... \n", request);
    send_status_req = zmq_send(requester, request, strlen(request), 0);

    // char received_data[max_floats * sizeof(float)];

    // Receive data via ZeroMQ
    char buffer[1024];  
    zmq_recv(requester, buffer, sizeof(buffer), 0);

    // Process received string
    int count = 0;
    const char delim[] = ";";
    char* token;
    float* float_array;

    // Count the number of floats
    char* copy = strdup(buffer);
    token = strtok(copy, delim);
    while (token != NULL) {
        count++;
        token = strtok(NULL, delim);
    }
    free(copy);

    // Allocate memory for the float array
    float_array = (float*)malloc(count * sizeof(float));

    // Convert string tokens to floats
    token = strtok(buffer, delim);
    for (int i = 0; i < count; ++i) {
        float_array[i] = strtof(token, NULL);
        token = strtok(NULL, delim);
    }

    // printf("Received float array:");
    // for (int i = 0; i < count; ++i) {
    //     printf(" %.4f", float_array[i]);
    // }
    // printf("\n");

    return float_array;
}




















    // const int max_data_size = 1024;

    // Receive data via ZeroMQ
    // char* received_data = (char*)malloc(max_data_size * sizeof(char));

    // Receive data via ZeroMQ
    // int recv_size = zmq_recv(requester, received_data, max_data_size - 1, 0);
    
    // if (recv_size > 0) {
    //     received_data[recv_size] = '\0'; // Null-terminate at the correct position
    //     printf("Received data: %s\n", received_data);

    //     return received_data;
    // } else {
    //     fprintf(stderr, "Error receiving data\n");
    //     free(received_data);
    //     return NULL;
    // }
    // }




//     float *received_value = (float *)malloc(sizeof(float));

//     // int recv_size = zmq_recv(requester, (char *)received_value, max_floats * sizeof(float), 0);
//     int recv_size = zmq_recv(requester, (char *)&received_value, sizeof(float *), 0);

//     if (recv_size > 0 && recv_size == sizeof(float *)) {
//         printf("Received pointer to array\n");
//         return received_value; // Return pointer to the received float array
//     } 
//     else 
//     {
//         printf("C: Error receiving float array\n");
//         return NULL;
//     }
// }

//     free(socket_address);
//     free(request);

//     return NULL; // Return NULL to indicate failure
// }
