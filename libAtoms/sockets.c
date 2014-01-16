#include <sys/socket.h>
#include <sys/types.h>
#include <netinet/in.h>
#include <netdb.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <errno.h>
#include <arpa/inet.h> 

#define MSG_LEN_SIZE 8
#define MSG_END_MARKER "done."
#define MSG_END_MARKER_SIZE strlen(MSG_END_MARKER)

int quip_recv_data(char *ip, int port, int client_id, char *data, int *data_len)
{
    int sockfd = 0, n = 0;
    char id_str[MSG_LEN_SIZE+1], msg_len_buff[MSG_LEN_SIZE+1], marker[MSG_END_MARKER_SIZE+1];
    int msg_len;
    int sent, totalsent, received, totalreceived;

    struct sockaddr_in serv_addr; 

    if((sockfd = socket(AF_INET, SOCK_STREAM, 0)) < 0)
    {
        printf("Could not create socket \n");
        return 1;
    } 

    memset(&serv_addr, 0, sizeof(serv_addr)); 

    serv_addr.sin_family = AF_INET;
    serv_addr.sin_port = htons(port);

    if(inet_pton(AF_INET, ip, &serv_addr.sin_addr)<=0)
    {
        printf("\n inet_pton error occured\n");
        return 1;
    } 

    if( connect(sockfd, (struct sockaddr *)&serv_addr, sizeof(serv_addr)) < 0)
    {
       printf("Connect Failed \n");
       return 1;
    }

    /* say hello and ask server for some Atoms (status 'A')
       identify ourselves with an ID number, formatted as an 8-byte ascii string */
    sprintf(id_str, "A%7d", client_id);

    totalsent = 0;
    while (totalsent < MSG_LEN_SIZE) {
      sent = send(sockfd, id_str+totalsent, MSG_LEN_SIZE - totalsent, 0);
      if (sent == 0) {
	printf("socket connection broken while sending client ID\n");
	return 1;
      }
      totalsent += sent;
    }

    /* now receive the length of the data, formatted as an 8-byte ascii string */
    memset(msg_len_buff, 0, sizeof(msg_len_buff));
    totalreceived = 0;
    while (totalreceived < MSG_LEN_SIZE)
    {
      received = recv(sockfd, msg_len_buff+totalreceived, MSG_LEN_SIZE-totalreceived, 0);
      if (received == 0) {
	printf("socket connection broken while reading length\n");
	return 1;
      }
      totalreceived += received;
    }
    sscanf(msg_len_buff, "%d", &msg_len);

    if (msg_len > *data_len) {
      printf("data to be sent is too large for receiver buffer\n");
      return 1;
    }
    *data_len = msg_len; // return the actual size of the data string

    /* now receive the data itself */
    memset(data, 0, sizeof(data));
    totalreceived = 0;
    while (totalreceived < msg_len)
    {
      received = recv(sockfd, data+totalreceived, msg_len-totalreceived, 0);
      if (received == 0) {
	printf("socket connection broken while reading data\n");
	return 1;
      }
      totalreceived += received;
    }

    /* and finally wait to receive the end marker */
    memset(marker, 0, sizeof(marker));
    totalreceived = 0;
    while (totalreceived < MSG_END_MARKER_SIZE)
    {
      received = recv(sockfd, marker+totalreceived, MSG_END_MARKER_SIZE-totalreceived, 0);
      if (received == 0) {
	printf("socket connection broken while reading data\n");
	return 1;
      }
      totalreceived += received;
    }
    
    close(sockfd);
    return 0;
}

int quip_send_data(char *ip, int port, int client_id, char *data, int data_len)
{
    int sockfd = 0, n = 0;
    char id_str[MSG_LEN_SIZE+1], msg_len_buff[MSG_LEN_SIZE+1], marker[MSG_END_MARKER_SIZE+1];
    int msg_len;
    int sent, totalsent, received, totalreceived;

    struct sockaddr_in serv_addr;

    if((sockfd = socket(AF_INET, SOCK_STREAM, 0)) < 0)
    {
        printf("Could not create socket \n");
        return 1;
    }

    memset(&serv_addr, 0, sizeof(serv_addr));

    serv_addr.sin_family = AF_INET;
    serv_addr.sin_port = htons(port);

    if(inet_pton(AF_INET, ip, &serv_addr.sin_addr)<=0)
    {
        printf("\n inet_pton error occured\n");
        return 1;
    }

    if( connect(sockfd, (struct sockaddr *)&serv_addr, sizeof(serv_addr)) < 0)
    {
       printf("Connect Failed \n");
       return 1;
    }

    /* say hello and tell server we've got some results (status 'R')
       identify ourselves with an ID number, formatted as an 8-byte ascii string */
    sprintf(id_str, "R%7d", client_id);

    totalsent = 0;
    while (totalsent < MSG_LEN_SIZE) {
      sent = send(sockfd, id_str+totalsent, MSG_LEN_SIZE - totalsent, 0);
      if (sent == 0) {
	printf("socket connection broken while sending client ID\n");
	return 1;
      }
      totalsent += sent;
    }

    /* Now send the length of the data, as an 8-byte string */
    sprintf(id_str, "%8d", data_len);

    totalsent = 0;
    while (totalsent < MSG_LEN_SIZE) {
      sent = send(sockfd, id_str+totalsent, MSG_LEN_SIZE - totalsent, 0);
      if (sent == 0) {
	printf("socket connection broken while sending data_len\n");
	return 1;
      }
      totalsent += sent;
    }

    /* send the data string itself */
    totalsent = 0;
    while (totalsent < data_len) {
      sent = send(sockfd, data+totalsent, data_len - totalsent, 0);
      if (sent == 0) {
	printf("socket connection broken while sending data\n");
	return 1;
      }
      totalsent += sent;
    }

    /* and finally wait to receive the end marker */
    memset(marker, 0, sizeof(marker));
    totalreceived = 0;
    while (totalreceived < MSG_END_MARKER_SIZE)
    {
      received = recv(sockfd, marker+totalreceived, MSG_END_MARKER_SIZE-totalreceived, 0);
      if (received == 0) {
	printf("socket connection broken while reading data\n");
	return 1;
      }
      totalreceived += received;
    }

    close(sockfd);
    return 0;
}
